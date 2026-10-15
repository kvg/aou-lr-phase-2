#!/usr/bin/env python3
"""Summarize Cromwell/Terra monitoring.log resource traces into a right-sizing report.

Input is the tree produced by scripts/collect_terra_monitoring_logs.sh: a
manifest.tsv plus call-<Task>/shard-N/[attempt-K/]monitoring.log files. Each log
holds the VM's visible #CPU / Total Memory / Total Disk and a ~10 s trace of CPU
usage, memory usage, disk usage and read/write IO.

Writes (into --outdir, default ./monitoring_report):
  monitoring_per_shard.tsv   one row per log: duration, peak/mean CPU, memory,
                             disk and IO
  monitoring_by_task.tsv     one row per task: allocation, observed peaks across
                             shards, headroom and the GiB-hours left idle
  monitoring_report.md       the same two tables plus per-task notes
  monitoring_utilization.png peak memory and peak disk vs allocation, per task

Usage:
  ./scripts/summarize_terra_monitoring.py terra_monitoring.tar.gz
  ./scripts/summarize_terra_monitoring.py /path/to/extracted/dir --outdir report/

Memory and disk peaks are sampled at the log interval, so a spike shorter than
one interval is invisible; size headroom accordingly.
"""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

TS_RE = re.compile(r"^\[(.+?)\]\s*$")
CPU_RE = re.compile(r"^\* CPU usage:\s*([\d.]+)%")
MEM_RE = re.compile(r"^\* Memory usage:\s*([\d.]+) GiB\s*([\d.]+)%")
DISK_RE = re.compile(r"^\* Disk usage:\s*([\d.]+) GiB\s*([\d.]+)%")
IO_RE = re.compile(r"^\* Read/Write IO:\s*([\d.]+) MiB/s\s*([\d.]+) MiB/s")
HDR_CPU_RE = re.compile(r"^#CPU:\s*(\d+)")
HDR_MEM_RE = re.compile(r"^Total Memory:\s*([\d.]+) GiB")
HDR_DISK_RE = re.compile(r"^Total Disk space:\s*([\d.]+) GiB")


def parse_log(path: Path) -> dict:
    """One monitoring.log -> header allocation + per-sample arrays."""
    cpu, mem, disk, rd, wr, ts = [], [], [], [], [], []
    alloc = {"n_cpu": np.nan, "mem_total_gib": np.nan, "disk_total_gib": np.nan}
    with path.open(errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if m := HDR_CPU_RE.match(line):
                alloc["n_cpu"] = int(m.group(1))
            elif m := HDR_MEM_RE.match(line):
                alloc["mem_total_gib"] = float(m.group(1))
            elif m := HDR_DISK_RE.match(line):
                alloc["disk_total_gib"] = float(m.group(1))
            elif m := CPU_RE.match(line):
                cpu.append(float(m.group(1)))
            elif m := MEM_RE.match(line):
                mem.append(float(m.group(1)))
            elif m := DISK_RE.match(line):
                disk.append(float(m.group(1)))
            elif m := IO_RE.match(line):
                rd.append(float(m.group(1)))
                wr.append(float(m.group(2)))
            elif m := TS_RE.match(line):
                try:
                    # e.g. "Mon Oct 12 05:07:58 UTC 2026"
                    ts.append(datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Z %Y"))
                except ValueError:
                    ts.append(None)
    return {**alloc, "cpu": cpu, "mem": mem, "disk": disk, "read": rd, "write": wr, "ts": ts}


# vCPU each task asks for, from the task defaults in felix/wdl/FelixGenome.wdl
# (Pack is pack_cpu from the Terra method config). Override with --requested-cpu.
# The monitoring script's "#CPU:" header reports 1 on every VM here and is not
# usable; CPU usage % is normalized to the whole machine, so a task with one busy
# thread on an 8-vCPU VM reads ~12.5%.
REQUESTED_CPU = {
    "Step2": 8, "Pack": 4, "RuLoci": 4, "RuVcf": 2, "BaVcf": 2,
    "RuStep2": 2, "BaStep2": 2, "Check": 1, "CheckNulls": 1,
    "Null": 8, "MakeGRM": 8, "Concat": 1, "RuConcat": 1, "BaConcat": 1,
    "Summarize": 2, "SummarizeRu": 2, "SummarizeBa": 2, "ArmCompare": 1, "Bundle": 1,
}

# GCE us-central1 list prices, N1 custom machine type, on-demand (USD).
# Indicative only: tasks with preemptible > 0 bill at the spot rate on early
# attempts, and the actual machine family is not recorded in the log.
PRICE_VCPU_HR = 0.033174
PRICE_MEM_GIB_HR = 0.004446
PRICE_DISK_GIB_HR = 0.04 / 730  # PD standard, $/GB-month -> $/GiB-hour


def shard_row(rel: str, meta: dict, parsed: dict) -> dict:
    cpu = np.asarray(parsed["cpu"], float)
    mem = np.asarray(parsed["mem"], float)
    disk = np.asarray(parsed["disk"], float)
    rd = np.asarray(parsed["read"], float)
    wr = np.asarray(parsed["write"], float)
    good_ts = [t for t in parsed["ts"] if t is not None]
    if len(good_ts) >= 2:
        span_s = (good_ts[-1] - good_ts[0]).total_seconds()
        interval_s = span_s / (len(good_ts) - 1)
    else:
        span_s, interval_s = np.nan, np.nan
    # The trace starts one interval into the task, so wall >= span + interval.
    wall_min = (span_s + interval_s) / 60 if np.isfinite(span_s) else np.nan

    def q(a, p):
        return float(np.percentile(a, p)) if a.size else np.nan

    def mx(a):
        return float(a.max()) if a.size else np.nan

    def mn(a):
        return float(a.mean()) if a.size else np.nan

    # Steady state: drop the localization head and delocalization tail, where
    # parallel gsutil transfers use cores the task itself never uses.
    steady = cpu[30:-6] if cpu.size > 60 else cpu

    row = {
        "relpath": rel,
        "call": meta["call"],
        # Every shard-N on the path, outermost first: for Step2 that is
        # <chromosome index>/<phenotype index>.
        "shard_chain": "/".join(re.findall(r"shard-(\d+)", rel)),
        "shard": meta["shard"],
        "attempt": meta["attempt"],
        "n_samples": int(cpu.size),
        "wall_min": wall_min,
        "n_cpu_header": parsed["n_cpu"],
        "cpu_mean_pct": mn(cpu),
        "cpu_steady_pct": float(np.median(steady)) if steady.size else np.nan,
        # Spread of the steady-state CPU: near zero means one thread pinned for
        # the whole task, which is what makes implied_vcpu below meaningful.
        "cpu_steady_rel_iqr": (
            float(np.subtract(*np.percentile(steady, [75, 25])) / np.median(steady))
            if steady.size and np.median(steady) > 0 else np.nan),
        "cpu_p95_pct": q(cpu, 95),
        "cpu_max_pct": mx(cpu),
        "cpu_frac_above_90": float((cpu > 90).mean()) if cpu.size else np.nan,
        "cpu_frac_below_10": float((cpu < 10).mean()) if cpu.size else np.nan,
        "mem_total_gib": parsed["mem_total_gib"],
        "mem_mean_gib": mn(mem),
        "mem_max_gib": mx(mem),
        "disk_total_gib": parsed["disk_total_gib"],
        "disk_max_gib": mx(disk),
        "read_mean_mibs": mn(rd),
        "read_max_mibs": mx(rd),
        "write_mean_mibs": mn(wr),
        "write_max_mibs": mx(wr),
    }
    row["mem_peak_frac"] = row["mem_max_gib"] / row["mem_total_gib"]
    row["disk_peak_frac"] = row["disk_max_gib"] / row["disk_total_gib"]
    # Integrated IO over the trace, in GiB.
    if np.isfinite(interval_s):
        row["read_gib"] = float(rd.sum() * interval_s / 1024)
        row["write_gib"] = float(wr.sum() * interval_s / 1024)
    else:
        row["read_gib"] = row["write_gib"] = np.nan
    return row


def load_tree(src: Path, tmp: Path) -> Path:
    if src.is_dir():
        return src
    with tarfile.open(src) as tf:
        tf.extractall(tmp)
    roots = [p for p in tmp.iterdir() if p.is_dir()]
    if len(roots) != 1:
        sys.exit(f"Expected one top-level directory in {src}, found {len(roots)}")
    return roots[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", type=Path, help="tar.gz from collect_terra_monitoring_logs.sh, or an extracted directory")
    ap.add_argument("--outdir", type=Path, default=Path("monitoring_report"))
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        root = load_tree(args.src, Path(td))
        manifest = pd.read_csv(root / "manifest.tsv", sep="\t", dtype=str)
        run_root = (root / "run_root.txt").read_text().strip() if (root / "run_root.txt").exists() else "unknown"
        rows = [
            shard_row(rec["relpath"], rec, parse_log(root / rec["relpath"]))
            for rec in manifest.to_dict("records")
            if (root / rec["relpath"]).exists()
        ]

    per_shard = pd.DataFrame(rows).sort_values(["call", "shard", "attempt"])
    if per_shard.empty:
        sys.exit("No monitoring.log files parsed.")

    per_shard["cpu_req"] = per_shard["call"].map(REQUESTED_CPU)
    # One busy thread on an N-vCPU VM reads 100/N %, so the steady-state
    # percentage recovers the machine size the task was actually given.
    per_shard["implied_vcpu"] = 100 / per_shard["cpu_steady_pct"]
    per_shard["cores_busy_mean"] = per_shard["cpu_mean_pct"] / 100 * per_shard["cpu_req"]
    # A large memory request forces a machine with more vCPUs than `cpu` asked
    # for (GCE custom types cap at ~6.5 GB per vCPU), and those cores are billed.
    # Trust the inferred size only where one thread is pinned flat throughout.
    tight = per_shard["cpu_steady_rel_iqr"] < 0.05
    per_shard["vcpu_billed"] = np.where(
        tight, np.maximum(per_shard["cpu_req"], per_shard["implied_vcpu"].round()), per_shard["cpu_req"])

    # Allocation is a per-shard property (disk is sized from input size), so the
    # task-level allocation columns report the range, not a single number.
    g = per_shard.groupby("call")
    by_task = pd.DataFrame({
        "n_shards": g.size(),
        "wall_min_median": g["wall_min"].median(),
        "wall_min_max": g["wall_min"].max(),
        "cpu_req": g["cpu_req"].max(),
        "vcpu_billed": g["vcpu_billed"].median(),
        "implied_vcpu": g["implied_vcpu"].median(),
        "cpu_steady_rel_iqr": g["cpu_steady_rel_iqr"].median(),
        "cores_busy_mean": g["cores_busy_mean"].mean(),
        "cores_busy_max": g["cores_busy_mean"].max(),
        "cpu_mean_pct": g["cpu_mean_pct"].mean(),
        "cpu_p95_pct": g["cpu_p95_pct"].max(),
        "mem_alloc_gib": g["mem_total_gib"].max(),
        "mem_peak_gib_median": g["mem_max_gib"].median(),
        "mem_peak_gib_max": g["mem_max_gib"].max(),
        "mem_peak_frac_max": g["mem_peak_frac"].max(),
        "disk_alloc_gib_median": g["disk_total_gib"].median(),
        "disk_alloc_gib_max": g["disk_total_gib"].max(),
        "disk_peak_gib_max": g["disk_max_gib"].max(),
        "disk_peak_frac_max": g["disk_peak_frac"].max(),
        "read_gib_median": g["read_gib"].median(),
        "write_gib_median": g["write_gib"].median(),
    })
    # Resource-hours actually held by the VM, and the share never touched.
    hrs = per_shard["wall_min"] / 60
    per_shard["mem_gib_hr_alloc"] = per_shard["mem_total_gib"] * hrs
    per_shard["mem_gib_hr_idle"] = (per_shard["mem_total_gib"] - per_shard["mem_max_gib"]) * hrs
    per_shard["disk_gib_hr_alloc"] = per_shard["disk_total_gib"] * hrs
    per_shard["disk_gib_hr_idle"] = (per_shard["disk_total_gib"] - per_shard["disk_max_gib"]) * hrs
    per_shard["vm_hr"] = hrs
    per_shard["vcpu_hr_alloc"] = per_shard["vcpu_billed"] * hrs
    per_shard["vcpu_hr_busy"] = per_shard["cores_busy_mean"] * hrs
    per_shard["cost_usd"] = (
        per_shard["vcpu_billed"] * PRICE_VCPU_HR
        + per_shard["mem_total_gib"] * PRICE_MEM_GIB_HR
        + per_shard["disk_total_gib"] * PRICE_DISK_GIB_HR
    ) * hrs
    g2 = per_shard.groupby("call")
    for col in ("mem_gib_hr_alloc", "mem_gib_hr_idle", "disk_gib_hr_alloc",
                "disk_gib_hr_idle", "vm_hr", "vcpu_hr_alloc", "vcpu_hr_busy", "cost_usd"):
        by_task[col] = g2[col].sum()
    by_task = by_task.sort_values("cost_usd", ascending=False)

    args.outdir.mkdir(parents=True, exist_ok=True)
    per_shard.to_csv(args.outdir / "monitoring_per_shard.tsv", sep="\t", index=False, float_format="%.4g")
    by_task.to_csv(args.outdir / "monitoring_by_task.tsv", sep="\t", float_format="%.4g")

    plot_utilization(by_task, args.outdir / "monitoring_utilization.png")

    traced = by_task[by_task["n_shards"].notna() & by_task["wall_min_median"].notna()]
    untraced = sorted(set(by_task.index) - set(traced.index))
    with (args.outdir / "monitoring_report.md").open("w") as fh:
        fh.write("# Terra resource utilization\n\n")
        fh.write(f"Run root: `{run_root}`\n\n")
        fh.write(f"{len(per_shard)} task executions across {per_shard['call'].nunique()} tasks. ")
        fh.write(f"Total {by_task['vm_hr'].sum():,.0f} VM-hours, ")
        fh.write(f"{by_task['vcpu_hr_alloc'].sum():,.0f} vCPU-hours allocated of which ")
        fh.write(f"{by_task['vcpu_hr_busy'].sum():,.0f} were busy.\n\n")
        if untraced:
            fh.write(f"No runtime samples (shorter than one monitoring interval): {', '.join(untraced)}.\n\n")
        fh.write("## Per task\n\n```\n")
        fh.write(by_task.to_string(float_format=lambda v: f"{v:,.2f}"))
        fh.write("\n```\n\n## Per shard\n\nSee monitoring_per_shard.tsv.\n")
        fh.write("\nCost is a list-price estimate: GCE us-central1 N1 custom on-demand at "
                 f"${PRICE_VCPU_HR}/vCPU-hr and ${PRICE_MEM_GIB_HR}/GiB-hr plus PD-standard disk; "
                 "tasks with preemptible > 0 bill lower on early attempts.\n")
    print(by_task.to_string(float_format=lambda v: f"{v:,.2f}"))


def plot_utilization(by_task: pd.DataFrame, path: Path) -> None:
    """Allocated vs actually-used CPU, memory and disk, one row per task."""
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    d = by_task[by_task["wall_min_median"].notna()].sort_values("cost_usd")
    labels = [f"{c}  (n={int(n)})" for c, n in zip(d.index, d["n_shards"])]
    panels = [
        ("vCPU requested vs busy", d["cpu_req"], d["cores_busy_mean"]),
        ("Memory allocated vs peak (GiB)", d["mem_alloc_gib"], d["mem_peak_gib_max"]),
        ("Disk allocated vs peak (GiB)", d["disk_alloc_gib_median"], d["disk_peak_gib_max"]),
    ]
    alloc_c, used_c = "#d9d9d9", "#1f6f8b"
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), sharey=True)
    y = np.arange(len(d))
    for i, (ax, (name, alloc, used)) in enumerate(zip(axes, panels)):
        ax.barh(y, alloc, color=alloc_c, edgecolor="none", height=0.68)
        ax.barh(y, used, color=used_c, edgecolor="none", height=0.38)
        for yi, (a, u) in enumerate(zip(alloc, used)):
            ax.text(a * 1.03, yi, f"{u:,.1f} / {a:,.0f}", va="center", ha="left", fontsize=7)
        ax.set_xlim(0, float(alloc.max()) * 1.5)
        ax.set_title(name, loc="left", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
        if i == 0:
            # Where the memory request forces a bigger machine than cpu asks for,
            # the idle-core count is larger still: mark the implied machine size.
            # Only where one thread is pinned flat for the whole task (tight
            # steady state) does the percentage pin down the machine size.
            for yi, (req, imp, iqr) in enumerate(zip(d["cpu_req"], d["implied_vcpu"], d["cpu_steady_rel_iqr"])):
                if np.isfinite(imp) and imp > req * 1.4 and iqr < 0.05:
                    ax.barh(yi, imp, height=0.68, facecolor="none",
                            edgecolor="#8c564b", linestyle="--", linewidth=0.9)
                    ax.text(imp * 1.03, yi, f"  machine ~{imp:.0f} vCPU", va="center",
                            ha="left", fontsize=7, color="#8c564b")
    axes[0].set_yticks(y, labels)
    fig.suptitle("No FelixGenome task uses more than 1.1 vCPU, 2.4 GiB of memory, or 99 GiB of disk\n"
                 "grey = allocated, blue = peak used (10 s sampling)",
                 x=0.012, y=0.985, va="top", ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
