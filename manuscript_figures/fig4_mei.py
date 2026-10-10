"""Figure 4: MEI catalog, SVA internal structure, and regulatory associations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.pyplot as plt

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    MEI_COLORS,
    OKABE,
    P2_ANC_N,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)


def _panel_mei_catalog(ax, rng: np.random.Generator) -> None:
    classes = ["Alu", "LINE-1", "SVA"]
    # Ancestry-stratified non-reference insertion counts (placeholder).
    # Totals ~ 12k Alu, 2.4k L1, 1.1k SVA polymorphic non-ref sites.
    base = {
        "Alu": 4200,
        "LINE-1": 780,
        "SVA": 360,
    }
    # Per-ancestry relative burden (AFR highest for L1/SVA).
    rel = {
        "AFR": (1.35, 1.55, 1.40),
        "AMR": (1.10, 1.05, 1.08),
        "EAS": (0.72, 0.60, 0.55),
        "EUR": (0.80, 0.70, 0.62),
        "MID": (0.90, 0.85, 0.78),
        "SAS": (0.88, 0.82, 0.90),
        "OTH": (1.05, 1.00, 1.02),
    }
    x = np.arange(len(ANC_ORDER) - 1)  # drop OTH for readability
    ancs = [a for a in ANC_ORDER if a != "OTH"]
    width = 0.26
    for i, cls in enumerate(classes):
        vals = [base[cls] * rel[a][i] * (P2_ANC_N[a] / 2000) ** 0.15 for a in ancs]
        ax.bar(
            x + (i - 1) * width,
            vals,
            width=width,
            color=MEI_COLORS[cls],
            label=cls,
            zorder=2,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(ancs)
    ax.set_ylabel("Non-reference MEI sites (placeholder)")
    ax.set_title("MEI catalog by ancestry", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right")


def _panel_sva_structure(ax, rng: np.random.Generator) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    # Schematic of an SVA element: hexamer, Alu-like, VNTR, SINE-R, polyA
    parts = [
        (0.4, 1.1, "#56B4E9", "(CCCTCT)n"),
        (1.5, 1.8, "#0072B2", "Alu-like"),
        (3.3, 2.6, "#E69F00", "VNTR"),
        (5.9, 2.4, "#009E73", "SINE-R"),
        (8.3, 1.1, "#D55E00", "polyA"),
    ]
    y = 4.55
    ax.text(0.4, 5.55, "Canonical SVA structure (not to scale)", fontsize=6.5, fontweight="medium")
    for x, w, color, lab in parts:
        ax.add_patch(Rectangle((x, y), w, 0.7, facecolor=color, edgecolor="none", alpha=0.85))
        ax.text(x + w / 2, y + 0.35, lab, ha="center", va="center", fontsize=5.8, color="white")
    ax.annotate(
        "",
        xy=(9.6, y + 0.35),
        xytext=(0.3, y + 0.35),
        arrowprops=dict(arrowstyle="-", color="#CCCCCC", lw=0.0),
    )
    ax.plot([0.4, 9.4], [y - 0.15, y - 0.15], color="#888888", lw=0.6)
    ax.text(0.4, y - 0.45, "5'", fontsize=6)
    ax.text(9.2, y - 0.45, "3'", fontsize=6)

    # VNTR length distributions by ancestry
    ax.text(0.4, 3.55, "SVA VNTR length (repeat copies)", fontsize=6.5)
    xs = np.linspace(5, 55, 200)
    for anc, mu, sd in [
        ("AFR", 28, 7.5),
        ("EUR", 21, 5.0),
        ("EAS", 18, 4.2),
        ("AMR", 24, 6.0),
    ]:
        ys = np.exp(-0.5 * ((xs - mu) / sd) ** 2)
        ys = ys / ys.max() * 1.35
        ax.plot(0.4 + (xs - 5) / 55 * 9.0, 0.45 + ys, color=ANC_COLORS[anc], lw=1.15, label=anc)
        ax.fill_between(
            0.4 + (xs - 5) / 55 * 9.0,
            0.45,
            0.45 + ys,
            color=ANC_COLORS[anc],
            alpha=0.10,
            linewidth=0,
        )
    ax.legend(loc="upper right", ncol=4, fontsize=5.8, frameon=False, bbox_to_anchor=(0.98, 0.58))
    ax.plot([0.4, 9.4], [0.45, 0.45], color="#222222", lw=0.6)
    for tick, lab in [(0.4, "5"), (5.0, "30"), (9.4, "55")]:
        ax.plot([tick, tick], [0.38, 0.45], color="#222222", lw=0.5)
        ax.text(tick, 0.18, lab, ha="center", fontsize=5.5, color="#444444")


def _panel_phewas(ax, rng: np.random.Generator) -> None:
    cats = [
        ("Endocrine", 18),
        ("Hematologic", 14),
        ("Immune", 22),
        ("Neuro", 16),
        ("GI", 15),
        ("Cardio", 20),
        ("Respiratory", 12),
        ("Renal", 11),
        ("Neoplasm", 17),
        ("Sense organs", 10),
    ]
    x0 = 0
    xticks = []
    xlabels = []
    colors = [
        OKABE["blue"],
        OKABE["orange"],
        OKABE["green"],
        OKABE["pink"],
        OKABE["vermillion"],
        OKABE["sky"],
        "#6B6B6B",
        "#8E44AD",
        "#2C7A7B",
        "#B08900",
    ]
    hit_labels = {
        "Hematologic": ("HBA SVA", "filled"),
        "Immune": ("HLA Alu", "filled"),
        "Neuro": ("SVA VNTR", "open"),
        "GI": ("UGT1A1 L1", "filled"),
    }
    for i, (cat, n) in enumerate(cats):
        x = x0 + np.linspace(0.2, n - 0.2, n * 14)
        y = np.clip(rng.exponential(0.55, x.size), 0, 4.0)
        if cat in hit_labels:
            lab, kind = hit_labels[cat]
            hit_x = x0 + n * 0.55
            y = np.maximum(y, 9.2 * np.exp(-0.5 * ((x - hit_x) / 0.35) ** 2))
            if kind == "open":
                ax.scatter(
                    [hit_x],
                    [9.4],
                    s=22,
                    facecolors="white",
                    edgecolors="#111111",
                    linewidths=0.8,
                    zorder=4,
                )
            else:
                ax.scatter([hit_x], [9.4], s=18, c="#111111", zorder=4, linewidths=0)
            ax.text(hit_x, 9.9, lab, ha="center", va="bottom", fontsize=5.4, color="#111111")
        ax.scatter(x, y, s=5, c=colors[i], alpha=0.75, linewidths=0, rasterized=True)
        xticks.append(x0 + n / 2)
        xlabels.append(cat)
        x0 += n
    ax.axhline(7.3, color="#D55E00", ls="--", lw=0.7)
    ax.set_xlim(0, x0)
    ax.set_ylim(0, 11.4)
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title("MEI and SVA-VNTR PheWAS (placeholder)", loc="left", pad=2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=6)
    despine(ax)
    ax.text(
        0.99,
        0.95,
        "Filled: insertion presence   ·   Open: SVA VNTR length",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
        color="#555555",
    )


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(21)
    ax_panel("fig4_mei_catalog", 2.45, lambda ax: _panel_mei_catalog(ax, rng), left=0.10, right=0.98, top=0.88, bottom=0.16)
    ax_panel("fig4_sva_structure", 2.55, lambda ax: _panel_sva_structure(ax, rng), left=0.10, right=0.98, top=0.90, bottom=0.14)
    ax_panel("fig4_mei_phewas", 2.65, lambda ax: _panel_phewas(ax, rng), left=0.10, right=0.98, top=0.88, bottom=0.16)
    fig = new_figure(7.35)
    gs = GridSpec(
        3,
        1,
        figure=fig,
        left=0.09,
        right=0.98,
        top=0.96,
        bottom=0.10,
        hspace=0.48,
        height_ratios=[1.0, 1.05, 1.15],
    )
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    _panel_mei_catalog(ax_a, rng)
    _panel_sva_structure(ax_b, rng)
    _panel_phewas(ax_c, rng)

    stamp_mockup(fig)
    return save(fig, "fig4_mei")
