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
RED = "#A33B32"

_PROGRESS_SEGMENTS = (
    ("done", "Done", GREEN),
    ("running", "Running", BLUE),
    ("failed", "Failed", RED),
    ("not_started", "Not started", GRAY),
)


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


def plot_cohort_progress(
    out_dir: Path,
    *,
    run_id: str,
    catalog: str,
    counts: dict[str, dict[str, int]],
    fully_done: int,
    keep_n: int,
    inflight_vcpu: int,
    e2_quota: int,
    n_active_jobs: int,
    updated_at: str,
) -> list[Path]:
    """16:9 ExpansionHunter progress figure for a slide.

    ``counts`` maps stage name (``minicram``, ``genotype``) to
    done / running / failed / not_started sample counts. The figure is left
    open so a notebook can display it.
    """
    plt = _plt()
    _style(plt)
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    stages = ("minicram", "genotype")
    y = [1, 0]
    left = [0, 0]
    total = max(int(keep_n), 1)
    for key, label, color in _PROGRESS_SEGMENTS:
        widths = [int((counts.get(stage) or {}).get(key) or 0) for stage in stages]
        ax.barh(y, widths, left=left, height=0.52, color=color, label=label)
        for yi, width, start in zip(y, widths, left):
            if width >= total * 0.06:
                ax.text(
                    start + width / 2,
                    yi,
                    f"{width:,}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=12,
                )
        left = [a + b for a, b in zip(left, widths)]
    ax.set_yticks(y)
    ax.set_yticklabels(["Minicram", "Genotype"])
    ax.set_xlim(0, total)
    ax.set_xlabel("Samples")
    ax.tick_params(axis="y", length=0)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.02), fontsize=11)
    fig.text(0.12, 0.93, f"ExpansionHunter    {run_id}", fontsize=20, ha="left", va="center")
    fig.text(
        0.12,
        0.875,
        f"{catalog}    ·    {keep_n:,} samples    ·    {updated_at}",
        fontsize=12,
        color=GRAY,
        ha="left",
        va="center",
    )
    fig.text(
        0.12,
        0.07,
        f"Fully genotyped    {fully_done:,}  /  {keep_n:,}",
        fontsize=18,
        ha="left",
        va="center",
    )
    fig.text(
        0.12,
        0.03,
        f"{n_active_jobs:,} jobs in flight    ·    {inflight_vcpu:,}  /  {e2_quota:,} E2 vCPU",
        fontsize=12,
        color=GRAY,
        ha="left",
        va="center",
    )
    fig.subplots_adjust(left=0.14, right=0.96, top=0.78, bottom=0.20)
    return save_fig(fig, out_dir, f"eh_progress_{run_id}")


# Validated (dataviz validate_palette.js, light surface): adjacent ΔE >= 21 normal and deutan.
RUN_DONE = "#2E7D4F"
RUN_RUNNING = "#5B9BE0"
RUN_FAILED = "#C8463A"
RUN_TODO = "#D5D9DF"
COST_SPENT = "#2F5D8A"
COST_INFLIGHT = "#5B9BE0"

_RUN_SEGMENTS = (
    ("done", "Done", RUN_DONE, None),
    ("running", "Running", RUN_RUNNING, None),
    ("failed", "Failed, will retry", RUN_FAILED, None),
    ("gave_up", "Gave up", RUN_FAILED, "////"),
    ("not_started", "Not started", RUN_TODO, None),
)


def _stack(ax, y, parts, total, *, fmt):
    """Horizontal stacked bar with 2px white gaps; labels only where they fit."""
    left = 0.0
    for value, color, hatch, label in parts:
        if value <= 0:
            continue
        ax.barh(
            y, value, left=left, height=0.42, color=color, hatch=hatch,
            edgecolor="white", linewidth=2, label=label,
        )
        if value >= total * 0.08:
            dark = color not in (RUN_TODO,)
            ax.text(left + value / 2, y, fmt(value), ha="center", va="center",
                    color="white" if dark else "#1F2937", fontsize=12)
        left += value
    return left


def plot_run_progress(
    out_dir: Path,
    *,
    run_id: str,
    catalog: str,
    counts: dict[str, int],
    total: int,
    spent_usd: float,
    committed_usd: float,
    projected_usd: float,
    budget_usd: float,
    actual_usd: float | None,
    status_line: str,
    updated_at: str,
) -> list[Path]:
    """16:9 slide: sample states on top, estimated dollars against the budget below."""
    plt = _plt()
    _style(plt)
    fig, (ax_s, ax_c) = plt.subplots(2, 1, figsize=(13.33, 7.5), gridspec_kw={"height_ratios": [1, 1]})

    parts = [(int(counts.get(k) or 0), c, h, lbl) for k, lbl, c, h in _RUN_SEGMENTS]
    _stack(ax_s, 0, parts, max(total, 1), fmt=lambda v: f"{int(v):,}")
    ax_s.set_xlim(0, max(total, 1))
    ax_s.set_yticks([])
    ax_s.set_xlabel("Samples")
    ax_s.spines["left"].set_visible(False)
    ax_s.legend(ncol=5, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=11, handlelength=1.2)
    ax_s.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))

    inflight = max(committed_usd - spent_usd, 0.0)
    rest = max(projected_usd - committed_usd, 0.0)
    xmax = max(budget_usd, projected_usd, committed_usd) * 1.05 or 1.0
    cost_parts = [
        (spent_usd, COST_SPENT, None, "Spent (est.)"),
        (inflight, COST_INFLIGHT, None, "In flight (est.)"),
        (rest, RUN_TODO, "....", f"Rest of cohort (projected total ${projected_usd:,.0f})"),
    ]
    end = _stack(ax_c, 0, cost_parts, xmax, fmt=lambda v: f"${v:,.0f}")
    if end < xmax * 0.7:
        ax_c.text(end + xmax * 0.01, 0, f"\\${spent_usd:,.2f} spent  ·  \\${projected_usd:,.0f} projected",
                  ha="left", va="center", fontsize=12, color="#1F2937")
    ax_c.axvline(budget_usd, color="#1F2937", lw=2)
    ax_c.text(budget_usd, 0.32, f" budget ${budget_usd:,.0f}", ha="left" if budget_usd < xmax * 0.85 else "right",
              va="center", fontsize=11, color="#1F2937")
    ax_c.set_xlim(0, xmax)
    ax_c.set_yticks([])
    ax_c.set_xlabel("USD")
    ax_c.spines["left"].set_visible(False)
    ax_c.legend(ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=11, handlelength=1.2)
    ax_c.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"${v:,.0f}"))

    done = int(counts.get("done") or 0)
    fig.text(0.06, 0.95, f"ExpansionHunter    {run_id}", fontsize=20, ha="left", va="center")
    fig.text(0.06, 0.905, f"{catalog}    ·    {done:,} / {total:,} genotyped ({done / max(total, 1):.1%})    ·    {updated_at}",
             fontsize=12, color=GRAY, ha="left", va="center")
    actual = "not available yet" if actual_usd is None else f"${actual_usd:,.2f}"
    if len(status_line) > 110:
        status_line = status_line[:107].rstrip() + "..."
    status_line = status_line.replace("$", "\\$")
    fig.text(0.06, 0.035, f"Billed so far: {actual}    ·    {status_line}", fontsize=11, color=GRAY, ha="left", va="center")
    fig.subplots_adjust(left=0.06, right=0.97, top=0.80, bottom=0.12, hspace=0.95)
    return save_fig(fig, out_dir, f"eh_run_{run_id}")


def main() -> None:
    here = Path(__file__).resolve().parents[2] / "expansion_hunter" / "batch" / "figures"
    paths = plot_smoke_711(here)
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
