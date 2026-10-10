"""Figure S2: SV callset construction, filtering, and novelty."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)


def _box(ax, x, y, w, h, text, color) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.03,rounding_size=0.07",
            facecolor=color,
            edgecolor="none",
            alpha=0.16,
        )
    )
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.03,rounding_size=0.07",
            facecolor="none",
            edgecolor=color,
            linewidth=0.85,
        )
    )
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=6.0, color="#222222")


def _panel_pipeline(ax) -> None:
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4.4)
    ax.axis("off")
    ax.set_title("Ensemble SV callset (placeholder schematic)", loc="left", pad=2)
    _box(ax, 0.15, 2.55, 1.7, 1.15, "PBSV\n(alignment)", OKABE["blue"])
    _box(ax, 2.05, 2.55, 1.7, 1.15, "Sniffles2\n(alignment)", OKABE["sky"])
    _box(ax, 3.95, 2.55, 1.7, 1.15, "PAV\n(assembly)", OKABE["green"])
    _box(ax, 6.15, 2.55, 2.15, 1.15, "Intra-sample\nTruvari collapse", OKABE["orange"])
    _box(ax, 8.55, 2.55, 1.6, 1.15, "XGBoost\nfilter", OKABE["pink"])
    _box(ax, 10.35, 2.55, 1.5, 1.15, "Stringent\n/ lenient", "#6B6B6B")
    _box(ax, 2.4, 0.45, 2.4, 1.15, "Inter-sample\nTruvari merge", OKABE["vermillion"])
    _box(ax, 5.3, 0.45, 2.5, 1.15, "Kanpig\ncohort regenotype", OKABE["blue"])
    _box(ax, 8.3, 0.45, 3.2, 1.15, "Phased SV + SNV/indel\nreference panel", OKABE["green"])
    for x0, x1, y in [(1.85, 2.05, 3.12), (3.75, 3.95, 3.12), (5.65, 6.15, 3.12), (8.30, 8.55, 3.12), (10.15, 10.35, 3.12)]:
        ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.8, mutation_scale=7))
    ax.plot([9.35, 9.35, 3.6], [2.55, 1.85, 1.85], color="#444444", lw=0.8)
    ax.annotate("", xy=(2.4, 1.02), xytext=(3.6, 1.85), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.8, mutation_scale=7))
    ax.annotate("", xy=(5.3, 1.02), xytext=(4.8, 1.02), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.8, mutation_scale=7))
    ax.annotate("", xy=(8.3, 1.02), xytext=(7.8, 1.02), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.8, mutation_scale=7))
    ax.text(
        0.15,
        0.12,
        "Calibrated on HPRC / HGSVC3 assemblies  ·  XGBoost trained per sample on Kanpig features",
        fontsize=5.8,
        color="#555555",
        ha="left",
    )


def _panel_roc(ax, rng: np.random.Generator) -> None:
    fpr = np.linspace(0, 0.55, 80)
    xgb = np.clip(1 - np.exp(-8.5 * fpr ** 0.62) + rng.normal(0, 0.008, fpr.size), 0, 1)
    xgb[0] = 0.02
    ax.plot(fpr, np.clip(xgb, 0, 1), color=OKABE["blue"], lw=1.6, label="XGBoost (per-sample)")
    # Heuristic operating points from the methods (placeholder positions near the text).
    ax.scatter([0.46], [0.95], s=28, c="#888888", zorder=4, label="Unfiltered")
    ax.scatter([0.20], [0.91], s=28, c=OKABE["orange"], zorder=4, label="Lenient (TPR proxy 0.9)")
    ax.scatter([0.11], [0.81], s=28, c=OKABE["green"], zorder=4, label="Stringent (TPR proxy 0.7)")
    ax.scatter([0.28], [0.78], s=22, marker="s", c=OKABE["pink"], zorder=4, label="≥2 callers")
    ax.plot([0, 0.55], [0, 0.55], ls=":", color="#CCCCCC", lw=0.7)
    ax.set_xlim(0, 0.55)
    ax.set_ylim(0.45, 1.02)
    ax.set_xlabel("False-positive rate (HPRC dipcall, placeholder)")
    ax.set_ylabel("True-positive rate")
    ax.set_title("Filtering ROC vs caller-support heuristics", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="lower right", fontsize=5.5)


def _panel_pr(ax) -> None:
    labs = ["Unfiltered\nmerge", "Lenient", "Stringent"]
    prec = [0.54, 0.80, 0.89]
    rec = [0.95, 0.91, 0.81]
    x = np.arange(len(labs))
    w = 0.36
    ax.bar(x - w / 2, prec, width=w, color=OKABE["blue"], label="Precision", zorder=2)
    ax.bar(x + w / 2, rec, width=w, color=OKABE["orange"], label="Recall", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Against HPRC dipcall (GRCh38)")
    ax.set_title("Cohort-level precision / recall", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right")
    for i, (p, r) in enumerate(zip(prec, rec)):
        ax.text(i - w / 2, p + 0.02, f"{p:.2f}", ha="center", fontsize=5.5, color=OKABE["blue"])
        ax.text(i + w / 2, r + 0.02, f"{r:.2f}", ha="center", fontsize=5.5, color=OKABE["orange"])


def _panel_novelty(ax) -> None:
    ancs = [a for a in ANC_ORDER if a != "OTH"]
    # Placeholder: AFR contributes the largest share of SVs absent from Phase 1 / HPRC.
    shared = np.array([0.52, 0.61, 0.74, 0.71, 0.68, 0.70])
    novel = 1.0 - shared
    x = np.arange(len(ancs))
    ax.bar(x, shared, color="#D0D0D0", label="Shared with Phase 1 or HPRC", zorder=2)
    ax.bar(x, novel, bottom=shared, color=[ANC_COLORS[a] for a in ancs], label="Novel in Phase 2 (placeholder)", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(ancs)
    ax.set_ylabel("Fraction of SVs ≥50 bp")
    ax.set_title("Novelty of the Phase 2 SV catalog by ancestry", loc="left", pad=2)
    ax.set_ylim(0, 1.08)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.5)
    for i, n in enumerate(novel):
        ax.text(i, shared[i] + n / 2, f"{n:.0%}", ha="center", va="center", fontsize=5.6, color="white")


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(8)
    ax_panel("figS2_pipeline", 2.15, _panel_pipeline, left=0.04, right=0.98, top=0.88, bottom=0.08)
    ax_panel("figS2_roc", 2.85, lambda ax: _panel_roc(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS2_pr", 2.85, _panel_pr, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS2_novelty", 2.55, _panel_novelty, left=0.08, right=0.98, top=0.88, bottom=0.16)
    fig = new_figure(7.55)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.95,
        bottom=0.07,
        wspace=0.30,
        hspace=0.48,
        height_ratios=[1.05, 1.05, 1.05],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])

    _panel_pipeline(ax_a)
    _panel_roc(ax_b, rng)
    _panel_pr(ax_c)
    _panel_novelty(ax_d)

    stamp_mockup(fig, extra="ROC and novelty fractions follow the methods narrative, not a locked callset")
    return save(fig, "figS2_callset")
