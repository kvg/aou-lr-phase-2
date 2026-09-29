"""Slide-ready matplotlib figures for EH Batch pilot / smoke briefings."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import cost

# Sample 1000000, 711-locus degenerate catalog. Batch RUNNING windows from
# statusEvents (not queue). Output sizes from gcloud storage ls -l.
SMOKE_711 = {
    "sample_id": "1000000",
    "catalog": "711-locus degenerate GRCh38",
    "minicram_job": "eh-mini-260929-051752",
    "genotype_job": "eh-gt-260929-053643",
    "minicram_running_s": 460.0,  # 05:19:27 → 05:27:07
    "genotype_running_s": 62.0,  # 05:38:17 → 05:39:19
    "minicram_queue_s": 90.0,  # SCHEDULED → RUNNING
    "genotype_queue_s": 85.0,
    "minicram_bytes": 70_740_829,
    "minicrai_bytes": 3_888,
    "eh_json_bytes": 574_432,
    "eh_vcf_bytes": 158_335,
    "machine": "e2-standard-4",
}


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required for figures (pip install matplotlib)") from exc
    return plt


def _style(plt) -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.frameon": False,
            "savefig.dpi": 180,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
        }
    )


BLUE = "#2F5D8A"
ORANGE = "#C46B2D"
GRAY = "#6B7280"
GREEN = "#3D6B4F"


def save_fig(fig, out_dir: Path, stem: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "pdf"):
        path = out_dir / f"{stem}.{ext}"
        fig.savefig(path)
        paths.append(path)
    return paths


def _p90_line(ax, values: list[float], *, color=ORANGE) -> float:
    q = cost.quantile(values, 0.90)
    ax.axvline(q, color=color, ls="--", lw=1.2, label=f"p90 = {q:.2g}")
    return q


def plot_smoke_711(out_dir: Path) -> list[Path]:
    """Figures from the proven one-sample 711-locus Batch run (no GCS)."""
    plt = _plt()
    _style(plt)
    s = SMOKE_711
    written: list[Path] = []

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    stages = ["Minicram", "Genotype"]
    queue = [s["minicram_queue_s"] / 60, s["genotype_queue_s"] / 60]
    run = [s["minicram_running_s"] / 60, s["genotype_running_s"] / 60]
    ax.bar(stages, queue, color=GRAY, label="Queue (scheduled → running)")
    ax.bar(stages, run, bottom=queue, color=BLUE, label="Running (billed VM time)")
    ax.set_ylabel("Minutes")
    ax.set_title("Wall time · sample 1000000 · 711-locus degenerate catalog")
    ax.legend(loc="upper right")
    fig.text(
        0,
        -0.08,
        "Source: Batch statusEvents · "
        f"{s['minicram_job']} · {s['genotype_job']} · 2026-09-29",
        color=GRAY,
        fontsize=8,
    )
    written += save_fig(fig, out_dir, "eh_smoke_711_walltime")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    labels = ["minicram.cram", "EH.json", "EH.vcf"]
    mib = [
        s["minicram_bytes"] / (1024**2),
        s["eh_json_bytes"] / (1024**2),
        s["eh_vcf_bytes"] / (1024**2),
    ]
    ax.barh(labels, mib, color=[BLUE, GREEN, GRAY])
    ax.set_xlabel("MiB")
    ax.set_title("GCS outputs · sample 1000000 · not the Nearline bytes fetched")
    for i, v in enumerate(mib):
        ax.text(v + 0.4, i, f"{v:.2f}" if v < 10 else f"{v:.1f}", va="center", fontsize=10)
    fig.text(
        0,
        -0.08,
        "Source: gs://aou-lr-phase2-resources/batchRuns/expansion_hunter/1000000/ · 2026-09-29",
        color=GRAY,
        fontsize=8,
    )
    written += save_fig(fig, out_dir, "eh_smoke_711_outputs")
    plt.close(fig)

    mini_usd = cost.vm_cost_usd(seconds=s["minicram_running_s"], machine=s["machine"])
    gt_usd = cost.vm_cost_usd(seconds=s["genotype_running_s"], machine=s["machine"])
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.bar(["Minicram VM", "Genotype VM"], [mini_usd, gt_usd], color=[BLUE, GREEN])
    ax.set_ylabel("USD (us-central1 list, e2-standard-4)")
    ax.set_title("Compute cost estimate · excludes Nearline retrieval")
    ax.bar_label(ax.containers[0], fmt="$%.3f")
    fig.text(
        0,
        -0.08,
        "Nearline $ is total_network_bytes × $0.01/GiB from data_transfer_stats.tsv (not plotted).",
        color=GRAY,
        fontsize=8,
    )
    written += save_fig(fig, out_dir, "eh_smoke_711_vm_cost")
    plt.close(fig)
    return written


def plot_pilot(
    *,
    costs: list[dict[str, Any]],
    resources: dict[str, list[dict[str, float]]],
    compute: dict[str, Any] | None,
    out_dir: Path,
    budget_usd: float = 0.0,
    n_keep: int = 0,
    planning_usd: float = 0.0,
) -> list[Path]:
    plt = _plt()
    _style(plt)
    written: list[Path] = []
    n = len(costs)
    caption = f"Source: EH Batch pilot · n={n} · p90 lines dashed · planning numbers not invoices"

    if costs:
        totals = [c["total_usd"] for c in costs]
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.hist(totals, bins=min(20, max(5, n // 3)), color=BLUE, edgecolor="white")
        _p90_line(ax, totals)
        ax.axvline(cost.quantile(totals, 0.50), color=GRAY, ls=":", lw=1.2, label=f"p50 = {cost.quantile(totals, 0.50):.3g}")
        ax.set_xlabel("Estimated USD per sample (VM + Nearline)")
        ax.set_ylabel("Samples")
        ax.set_title("Pilot cost per sample")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, "eh_pilot_cost_hist")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        idx = list(range(n))
        mini = [c.get("minicram_vm_usd", 0) for c in costs]
        gt = [c.get("genotype_vm_usd", 0) for c in costs]
        near = [c.get("nearline_usd", 0) for c in costs]
        ax.bar(idx, mini, color=BLUE, label="Minicram VM")
        ax.bar(idx, gt, bottom=mini, color=GREEN, label="Genotype VM")
        bottom2 = [a + b for a, b in zip(mini, gt)]
        ax.bar(idx, near, bottom=bottom2, color=ORANGE, label="Nearline retrieval")
        ax.set_xlabel("Pilot sample (arbitrary order)")
        ax.set_ylabel("USD")
        ax.set_title("Cost stack per sample")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, "eh_pilot_cost_stack")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.scatter(
            [c["minicram_seconds"] / 60 for c in costs],
            [c["genotype_seconds"] / 60 for c in costs],
            c=BLUE,
            alpha=0.75,
            edgecolors="none",
        )
        ax.set_xlabel("Minicram running (min)")
        ax.set_ylabel("Genotype running (min)")
        ax.set_title("Stage wall time per sample")
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, "eh_pilot_walltime")
        plt.close(fig)

        gib = [c["network_bytes"] / (1024**3) for c in costs]
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.hist(gib, bins=min(20, max(5, n // 3)), color=ORANGE, edgecolor="white")
        _p90_line(ax, gib)
        ax.set_xlabel("Nearline bytes fetched (GiB)")
        ax.set_ylabel("Samples")
        ax.set_title("CRAM container retrieval (not minicram output size)")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, "eh_pilot_nearline")
        plt.close(fig)

    for stage, color in (("minicram", BLUE), ("genotype", GREEN)):
        rows = resources.get(stage) or []
        if len(rows) < 2:
            continue
        rss = [r["peak_rss_bytes"] / (1024**2) for r in rows]
        cores = [r["cpu_cores_avg"] for r in rows]
        rec = (compute or {}).get(stage) or {}
        req_mem = rec.get("memoryMib") or cost.DEFAULT_COMPUTE[stage]["memoryMib"]
        req_cpu = (rec.get("cpuMilli") or cost.DEFAULT_COMPUTE[stage]["cpuMilli"]) / 1000

        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.hist(rss, bins=min(16, max(5, len(rss) // 3)), color=color, edgecolor="white")
        _p90_line(ax, rss)
        ax.axvline(req_mem, color=GRAY, ls=":", label=f"request {req_mem} MiB")
        ax.set_xlabel("Peak RSS (MiB), process tree")
        ax.set_ylabel("Samples")
        ax.set_title(f"Peak memory · {stage}")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, f"eh_pilot_rss_{stage}")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        ax.hist(cores, bins=min(16, max(5, len(cores) // 3)), color=color, edgecolor="white")
        _p90_line(ax, cores)
        ax.axvline(req_cpu, color=GRAY, ls=":", label=f"request {req_cpu:.0f} vCPU")
        ax.set_xlabel("Average CPU cores (cpu-sec / wall-sec)")
        ax.set_ylabel("Samples")
        ax.set_title(f"CPU use · {stage}")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, f"eh_pilot_cpu_{stage}")
        plt.close(fig)

    if n_keep and planning_usd and budget_usd >= 0:
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        xs = list(range(0, n_keep + max(1, n_keep // 10), max(1, n_keep // 20))) or [n_keep]
        if xs[-1] != n_keep:
            xs.append(n_keep)
        ax.plot(xs, [x * planning_usd for x in xs], color=BLUE, label="n × p90 unit")
        if budget_usd > 0:
            ax.axhline(budget_usd, color=ORANGE, ls="--", label=f"budget ${budget_usd:,.0f}")
        ax.scatter([n_keep], [n_keep * planning_usd], color=BLUE, zorder=3)
        ax.set_xlabel("Samples on keep list")
        ax.set_ylabel("Estimated USD")
        ax.set_title("Budget vs keep-list size (p90 unit, no retry margin on this axis)")
        ax.legend()
        fig.text(0, -0.08, caption, color=GRAY, fontsize=8)
        written += save_fig(fig, out_dir, "eh_pilot_budget")
        plt.close(fig)

    return written


def main() -> None:
    here = Path(__file__).resolve().parents[2] / "expansion_hunter" / "batch" / "figures"
    paths = plot_smoke_711(here)
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
