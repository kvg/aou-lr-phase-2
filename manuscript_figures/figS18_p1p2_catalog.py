"""Figure S18: Phase 1 vs Phase 2 variant catalogs from Table 2.

Observed site counts (tab:callset). Not overlap/Truvari — that remains Fig. S6.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from .style import (
    OKABE,
    P1_N,
    P2_N,
    despine,
    format_millions,
    label_panel,
    new_figure,
    save,
    stamp_mockup,
)

C1 = "#B0B0B0"
C2 = OKABE["blue"]
C50 = OKABE["blue"]
C20 = OKABE["orange"]

# Table 2 (tab:callset).
P1_SNV, P2_SNV = 54_790_068, 245_289_070
P1_INS20, P2_INS20_SMALL = 3_932_966, 13_218_208  # <20 bp
P1_DEL20, P2_DEL20_SMALL = 4_734_797, 17_707_931
P1_DEL50, P2_DEL50 = 165_717, 599_220
P1_DEL_GE20, P2_DEL_GE20 = 165_724, 1_575_018
P1_INS50, P2_INS50 = 498_090, 1_875_567
P1_INS_GE20, P2_INS_GE20 = 498_094, 3_501_898
P1_TR, P2_TR = 1_784_804, 4_439_813
P2_BND, P2_LARGE = 130_740, 225_241

P1_DEL_2049 = P1_DEL_GE20 - P1_DEL50  # 7
P2_DEL_2049 = P2_DEL_GE20 - P2_DEL50
P1_INS_2049 = P1_INS_GE20 - P1_INS50  # 4
P2_INS_2049 = P2_INS_GE20 - P2_INS50

SAMPLE_FOLD = P2_N / P1_N

# ≥50 bp sequence-context percentages (US / RM / SD / SR / CMRG). CMRG is 0.
CTX = ["US", "RM", "SD", "SR"]
DEL_CTX_P1 = np.array([3.9, 58.6, 12.6, 82.4])
DEL_CTX_P2 = np.array([7.1, 66.6, 12.6, 69.5])
INS_CTX_P1 = np.array([3.6, 71.1, 8.3, 84.6])
INS_CTX_P2 = np.array([3.9, 66.0, 10.4, 84.7])


def _fold(a: float, b: float) -> float:
    return b / a


def _panel_catalog(ax) -> None:
    labels = [
        "SNVs",
        "INS <20 bp",
        "DEL <20 bp",
        r"DEL $\geq$50 bp",
        r"INS $\geq$50 bp",
        "TR loci",
        "Breakends",
        r"Large $>$10 kb",
    ]
    p1 = np.array([P1_SNV, P1_INS20, P1_DEL20, P1_DEL50, P1_INS50, P1_TR, 0.0, 1_997], dtype=float)
    p2 = np.array(
        [P2_SNV, P2_INS20_SMALL, P2_DEL20_SMALL, P2_DEL50, P2_INS50, P2_TR, P2_BND, P2_LARGE],
        dtype=float,
    )
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, np.maximum(p1, 0.8), width=w, color=C1, zorder=2, label=f"Phase 1  n={P1_N:,}")
    ax.bar(x + w / 2, p2, width=w, color=C2, zorder=2, label=f"Phase 2  n={P2_N:,}")
    ax.set_yscale("log")
    ax.set_ylim(0.7, 6e8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Catalog sites")
    ax.set_title("Table 2 catalogs  ·  unique sites, not per-genome medians", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.8)
    for i, (a, b) in enumerate(zip(p1, p2)):
        ax.text(i - w / 2, max(a, 1) * 1.12, format_millions(a) if a else "—", ha="center", va="bottom", fontsize=5.2, color="#555555")
        ax.text(i + w / 2, b * 1.12, format_millions(b), ha="center", va="bottom", fontsize=5.2, color="#222222")
    ax.text(
        0.01,
        0.06,
        "Breakends and >10 kb were not reported for Phase 1. TR catalogs differ (Adotto vs TRExplorer).",
        transform=ax.transAxes,
        fontsize=5.3,
        color="#666666",
    )


def _panel_2049(ax) -> None:
    """Phase 1 main callset was effectively ≥50 bp; Phase 2 adds a 20–49 bp slice."""
    labs = ["DEL", "INS"]
    p1_50 = np.array([P1_DEL50, P1_INS50], dtype=float)
    p1_mid = np.array([P1_DEL_2049, P1_INS_2049], dtype=float)
    p2_50 = np.array([P2_DEL50, P2_INS50], dtype=float)
    p2_mid = np.array([P2_DEL_2049, P2_INS_2049], dtype=float)
    x = np.arange(len(labs))
    w = 0.32
    ax.bar(x - w / 2, p1_50, width=w, color=C50, zorder=2)
    ax.bar(x - w / 2, p1_mid, width=w, bottom=p1_50, color=C20, zorder=3)
    ax.bar(x + w / 2, p2_50, width=w, color=C50, zorder=2)
    ax.bar(x + w / 2, p2_mid, width=w, bottom=p2_50, color=C20, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=8.0)
    ax.set_xlim(-0.7, 1.7)
    ax.set_ylabel("Sites in the main SV callset")
    ax.set_title("20–49 bp is new in Phase 2", loc="left", pad=2)
    despine(ax)
    ax.set_ylim(0, 4.05e6)
    for i, (lo, mid, phase, dx) in enumerate(
        [
            (p1_50[0], p1_mid[0], "P1", -w / 2),
            (p2_50[0], p2_mid[0], "P2", w / 2),
            (p1_50[1], p1_mid[1], "P1", -w / 2),
            (p2_50[1], p2_mid[1], "P2", w / 2),
        ]
    ):
        xi = (0 if i < 2 else 1) + dx
        ax.text(xi, 80_000, phase, ha="center", va="bottom", fontsize=6.0, color="white", fontweight="bold")
        if mid > 50_000:
            ax.text(xi, lo + mid / 2, format_millions(mid), ha="center", va="center", fontsize=5.6, color="#222222")
        ax.text(xi, lo + mid + 6e4, format_millions(lo + mid), ha="center", va="bottom", fontsize=5.8, color="#222222")
    ax.plot([], [], color=C50, lw=6, label=r"$\geq$50 bp")
    ax.plot([], [], color=C20, lw=6, label="20–49 bp")
    ax.legend(loc="upper left", fontsize=5.8)
    ax.text(
        0.98,
        0.08,
        "Phase 1 DEL 20–49 bp  n=7\nPhase 1 INS 20–49 bp  n=4",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.4,
        color="#555555",
    )


def _panel_fold(ax) -> None:
    rows = [
        ("TR loci*", _fold(P1_TR, P2_TR)),
        (r"DEL $\geq$50 bp", _fold(P1_DEL50, P2_DEL50)),
        ("SNVs", _fold(P1_SNV, P2_SNV)),
        (r"INS $\geq$50 bp", _fold(P1_INS50, P2_INS50)),
        ("DEL <20 bp", _fold(P1_DEL20, P2_DEL20_SMALL)),
        ("INS <20 bp", _fold(P1_INS20, P2_INS20_SMALL)),
        (r"INS $\geq$20 bp", _fold(P1_INS_GE20, P2_INS_GE20)),
        (r"DEL $\geq$20 bp", _fold(P1_DEL_GE20, P2_DEL_GE20)),
        ("Samples", SAMPLE_FOLD),
    ]
    y = np.arange(len(rows))
    folds = np.array([f for _, f in rows])
    colors = [OKABE["grey"] if lab == "Samples" else C2 for lab, _ in rows]
    ax.axvline(SAMPLE_FOLD, color=C1, lw=0.8, ls="--", zorder=1)
    ax.hlines(y, 0, folds, color=colors, lw=1.15, zorder=2)
    ax.scatter(folds, y, s=22, c=colors, zorder=3, edgecolors="white", linewidths=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    ax.set_xlabel("Phase 2 / Phase 1")
    ax.set_xlim(0, 13.2)
    ax.set_title("Catalog fold vs sample-size fold", loc="left", pad=2)
    despine(ax)
    for yi, f in zip(y, folds):
        ax.text(f + 0.18, yi, f"{f:.1f}x", va="center", fontsize=5.6, color="#333333")
    ax.text(
        SAMPLE_FOLD,
        len(rows) - 0.35,
        "n",
        ha="center",
        va="bottom",
        fontsize=5.4,
        color="#666666",
    )
    ax.text(0.02, 0.04, "*different TR catalogs", transform=ax.transAxes, fontsize=5.2, color="#666666")


def _panel_context(ax, p1: np.ndarray, p2: np.ndarray, title: str) -> None:
    x = np.arange(len(CTX))
    w = 0.36
    ax.bar(x - w / 2, p1, width=w, color=C1, zorder=2, label="Phase 1")
    ax.bar(x + w / 2, p2, width=w, color=C2, zorder=2, label="Phase 2")
    ax.set_xticks(x)
    ax.set_xticklabels(CTX)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of " + r"$\geq$50 bp sites")
    ax.set_title(title, loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.5)
    for i, (a, b) in enumerate(zip(p1, p2)):
        ax.text(i - w / 2, a + 1.5, f"{a:.1f}", ha="center", va="bottom", fontsize=5.2, color="#555555")
        ax.text(i + w / 2, b + 1.5, f"{b:.1f}", ha="center", va="bottom", fontsize=5.2, color="#222222")


def render() -> tuple[Path, Path]:
    fig = new_figure(8.35)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.09,
        right=0.985,
        top=0.96,
        bottom=0.055,
        wspace=0.32,
        hspace=0.42,
        height_ratios=[1.18, 1.12, 1.0],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_e = fig.add_subplot(gs[2, 1])

    _panel_catalog(ax_a)
    label_panel(ax_a, "A", x=-0.05, y=1.06)
    _panel_2049(ax_b)
    label_panel(ax_b, "B", x=-0.12, y=1.06)
    _panel_fold(ax_c)
    label_panel(ax_c, "C", x=-0.22, y=1.06)
    _panel_context(ax_d, DEL_CTX_P1, DEL_CTX_P2, r"Deletions $\geq$50 bp, sequence context")
    label_panel(ax_d, "D", x=-0.12, y=1.08)
    _panel_context(ax_e, INS_CTX_P1, INS_CTX_P2, r"Insertions $\geq$50 bp, sequence context")
    label_panel(ax_e, "E", x=-0.12, y=1.08)

    stamp_mockup(fig, extra="Counts are Table 2 observed site totals")
    return save(fig, "figS18_p1p2_catalog")
