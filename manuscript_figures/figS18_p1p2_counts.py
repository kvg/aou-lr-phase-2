"""Figure S18: Phase 1 vs Phase 2 catalog counts from Table 2.

Not a Truvari overlap (that is Fig. S6). This plate is the Table 2 numbers:
SNV/indel site counts, SV ≥50 vs 20–49 bp, sequence-context mix, and fold
change versus the 12-fold jump in n.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch

from .style import (
    OKABE,
    P1_N,
    P2_N,
    PANELDIR,
    despine,
    format_millions,
    label_panel,
    new_figure,
    save,
    stamp_mockup,
)

# Table 2 (tab:callset). Phase 1 mid-pass; Phase 2 combined PacBio discovery.
P1 = {
    "snv": 54_790_068,
    "ins_lt20": 3_932_966,
    "del_lt20": 4_734_797,
    "del_50": 165_717,
    "del_20": 165_724,
    "ins_50": 498_090,
    "ins_20": 498_094,
    "inv_50": 0,
    "large": 1_997,
    "tr": 1_784_804,
}
P2 = {
    "snv": 245_289_070,
    "ins_lt20": 13_218_208,
    "del_lt20": 17_707_931,
    "del_50": 599_220,
    "del_20": 1_575_018,
    "ins_50": 1_875_567,
    "ins_20": 3_501_898,
    "inv_50": 41,
    "bnd": 130_740,
    "large": 225_241,
    "tr": 4_439_813,
}

# ≥50 bp GIAB-style context percentages (Table 2). US / RM / SD / SR / CMRG.
CTX_DEL_P1 = np.array([3.9, 58.6, 12.6, 82.4, 0.0])
CTX_DEL_P2 = np.array([7.1, 66.6, 12.6, 69.5, 0.0])
CTX_INS_P1 = np.array([3.6, 71.1, 8.3, 84.6, 0.0])
CTX_INS_P2 = np.array([3.9, 66.0, 10.4, 84.7, 0.0])
CTX_LABS = ["US", "RM", "SD", "SR", "CMRG"]

C_P1 = "#B0B0B0"
C_P2 = OKABE["blue"]
C_SMALL = OKABE["orange"]  # 20–49 bp slice


def _panel_scale(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(
        f"Catalog size  ·  Phase 1 n={P1_N:,} AFR mid-pass  vs  Phase 2 n={P2_N:,} combined",
        loc="left",
        pad=2,
    )
    items = [
        (f"{P2['snv'] / P1['snv']:.1f}x", "SNVs"),
        (f"{P2['del_lt20'] / P1['del_lt20']:.1f}x", "DEL <20 bp"),
        (f"{P2['ins_lt20'] / P1['ins_lt20']:.1f}x", "INS <20 bp"),
        (f"{P2['del_50'] / P1['del_50']:.1f}x", "DEL >=50 bp"),
        (f"{P2['ins_50'] / P1['ins_50']:.1f}x", "INS >=50 bp"),
        (f"{P2_N / P1_N:.1f}x", "Sample n"),
    ]
    n = len(items)
    gap = 0.012
    w = (1.0 - gap * (n - 1)) / n
    for i, (num, lab) in enumerate(items):
        x0 = i * (w + gap)
        ax.add_patch(
            FancyBboxPatch(
                (x0, 0.08),
                w,
                0.84,
                boxstyle="round,pad=0.012,rounding_size=0.04",
                facecolor="#F4F4F4",
                edgecolor="#DDDDDD",
                linewidth=0.6,
                transform=ax.transAxes,
                clip_on=False,
            )
        )
        ax.text(x0 + w / 2, 0.62, num, ha="center", va="center", fontsize=12.0, fontweight="bold", transform=ax.transAxes)
        ax.text(x0 + w / 2, 0.28, lab, ha="center", va="center", fontsize=5.8, color="#444444", transform=ax.transAxes)


def _panel_counts(ax, *, title: bool = True) -> None:
    labels = [
        "SNVs",
        "DEL\n<20 bp",
        "INS\n<20 bp",
        "DEL\n$\\geq$50 bp",
        "INS\n$\\geq$50 bp",
        "TRGT\nloci",
    ]
    p1 = np.array([P1["snv"], P1["del_lt20"], P1["ins_lt20"], P1["del_50"], P1["ins_50"], P1["tr"]], dtype=float)
    p2 = np.array([P2["snv"], P2["del_lt20"], P2["ins_lt20"], P2["del_50"], P2["ins_50"], P2["tr"]], dtype=float)
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, p1 / 1e6, width=w, color=C_P1, zorder=2, label=f"Phase 1  n={P1_N:,}")
    ax.bar(x + w / 2, p2 / 1e6, width=w, color=C_P2, zorder=2, label=f"Phase 2  n={P2_N:,}")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.0, linespacing=1.05)
    ax.tick_params(axis="x", pad=2)
    ax.set_yscale("log")
    ax.set_ylim(0.08, 900)
    ax.set_ylabel("Catalog sites (millions)")
    if title:
        ax.set_title("Joint-callset site counts (Table 2)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.6)
    for i, (a, b) in enumerate(zip(p1, p2)):
        ax.text(i - w / 2, a / 1e6 * 1.12, format_millions(a), ha="center", va="bottom", fontsize=5.0, color="#555555", rotation=90)
        ax.text(i + w / 2, b / 1e6 * 1.12, format_millions(b), ha="center", va="bottom", fontsize=5.0, color=C_P2, rotation=90)


def _panel_svlen(ax) -> None:
    """Phase 1 ≥20 and ≥50 are almost the same; Phase 2 adds a real 20–49 bp class."""
    cats = ["DEL", "INS"]
    p1_50 = np.array([P1["del_50"], P1["ins_50"]], dtype=float)
    p1_2049 = np.array([P1["del_20"] - P1["del_50"], P1["ins_20"] - P1["ins_50"]], dtype=float)
    p2_50 = np.array([P2["del_50"], P2["ins_50"]], dtype=float)
    p2_2049 = np.array([P2["del_20"] - P2["del_50"], P2["ins_20"] - P2["ins_50"]], dtype=float)
    x = np.arange(len(cats))
    w = 0.32
    ax.bar(x - w / 2, p1_50 / 1e3, width=w, color=C_P1, zorder=2, label="Phase 1 ≥50 bp")
    ax.bar(x - w / 2, p1_2049 / 1e3, width=w, bottom=p1_50 / 1e3, color="#E2E2E2", zorder=2, label="Phase 1 20–49 bp")
    ax.bar(x + w / 2, p2_50 / 1e3, width=w, color=C_P2, zorder=2, label="Phase 2 ≥50 bp")
    ax.bar(x + w / 2, p2_2049 / 1e3, width=w, bottom=p2_50 / 1e3, color=C_SMALL, zorder=2, label="Phase 2 20–49 bp")
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("Sites (thousands)")
    ax.set_title("The 20–49 bp class is new in Phase 2", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.3)
    ax.text(
        0 - w / 2,
        (p1_50[0] + p1_2049[0]) / 1e3 + 40,
        f"20–49 bp  n={int(p1_2049[0])}",
        ha="center",
        va="bottom",
        fontsize=5.3,
        color="#666666",
    )
    ax.text(
        1 + w / 2,
        (p2_50[1] + p2_2049[1]) / 1e3 * 0.55,
        format_millions(p2_2049[1]),
        ha="center",
        va="center",
        fontsize=6.0,
        color="#222222",
        fontweight="bold",
    )
    ax.set_ylim(0, 4.05e3)


def _panel_context(ax) -> None:
    x = np.arange(len(CTX_LABS))
    w = 0.18
    ax.bar(x - 1.5 * w, CTX_DEL_P1, width=w, color="#D0D0D0", zorder=2, label="P1 DEL")
    ax.bar(x - 0.5 * w, CTX_DEL_P2, width=w, color=OKABE["blue"], zorder=2, label="P2 DEL")
    ax.bar(x + 0.5 * w, CTX_INS_P1, width=w, color="#E8C9A0", zorder=2, label="P1 INS")
    ax.bar(x + 1.5 * w, CTX_INS_P2, width=w, color=OKABE["orange"], zorder=2, label="P2 INS")
    ax.set_xticks(x)
    ax.set_xticklabels(CTX_LABS)
    ax.set_ylabel("% of ≥50 bp sites")
    ax.set_ylim(0, 100)
    ax.set_title("Sequence context (GIAB strata, ≥50 bp)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", ncol=2, fontsize=5.3)
    ax.text(0.02, 0.92, "RM, SD, SR overlap; not a partition", transform=ax.transAxes, fontsize=5.3, color="#666666")


def _panel_extra(ax) -> None:
    labels = ["Inversions\n≥50 bp", "Breakends", "Large events\n>10 kb", "TRGT loci"]
    p1 = np.array([P1["inv_50"], 0.0, P1["large"], P1["tr"]], dtype=float)
    p2 = np.array([P2["inv_50"], P2["bnd"], P2["large"], P2["tr"]], dtype=float)
    x = np.arange(len(labels))
    w = 0.36
    ax.bar(x - w / 2, np.where(p1 > 0, p1, np.nan), width=w, color=C_P1, zorder=2)
    ax.bar(x + w / 2, p2, width=w, color=C_P2, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.0)
    ax.set_yscale("log")
    ax.set_ylim(0.8, 1.2e7)
    ax.set_ylabel("Count")
    ax.set_title("Inversions, breakends, large events and TRGT", loc="left", pad=2)
    despine(ax)
    for i, (a, b) in enumerate(zip(p1, p2)):
        if a > 0:
            ax.text(i - w / 2, a * 1.15, format_millions(a), ha="center", va="bottom", fontsize=5.4, color="#555555")
        else:
            ax.text(i - w / 2, 1.2, "—", ha="center", va="bottom", fontsize=6.5, color="#888888")
        ax.text(i + w / 2, b * 1.12, format_millions(b) if b >= 1000 else f"{int(b)}", ha="center", va="bottom", fontsize=5.4, color=C_P2)


def _panel_fold_by_length(fig, cell) -> plt.Axes:
    """Phase 2 / Phase 1 sites per Fig. 2 length bin, ≤100 bp and 100 bp–10 kb."""
    from . import fig2_options as f2

    data = f2.spectrum_data(np.random.default_rng(0), region_path=f2._region_bins_path())
    facets = [f for f in data["facets"] if f["spec"]["layer"] == "integrated"]
    sub = cell.subgridspec(1, len(facets), width_ratios=[f["spec"]["width"] for f in facets], wspace=0.06)
    axes = [fig.add_subplot(sub[0, i]) for i in range(len(facets))]
    for i, (ax, f) in enumerate(zip(axes, facets)):
        spec = f["spec"]
        for d in ("INS", "DEL"):
            p2v, p1v = (f["ins"], f.get("p1_ins")) if d == "INS" else (f["del"], f.get("p1_del"))
            if p1v is None:
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                y = np.where(p1v > 0, p2v / p1v, np.nan)
            ax.step(f["edges"], np.r_[y, y[-1]], where="post", color=f2.FOLD_COLORS[d], lw=0.6)
            if i == 1:
                ax.text(2_500, 19 if d == "DEL" else 2.1, d, ha="center", va="center", fontsize=6.0, color=f2.FOLD_COLORS[d], fontweight="bold")
        ax.axhline(P2_N / P1_N, color="#777777", lw=0.5, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(1, 30)
        ax.set_yticks([1, 3, 10, 30])
        ax.set_yticklabels(["1", "3", "10", "30"])
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.grid(axis="y", color="#EEEEEE", lw=0.4)
        ax.set_xscale("log" if spec["xscale"] == "log" else "linear")
        ax.set_xlim(*spec["xlim"])
        ax.set_xticks(spec["xticks"])
        ax.set_xticklabels(spec.get("labels") or [f"{t:g}" for t in spec["xticks"]])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_title(spec["title"], fontsize=6.5, pad=2)
        despine(ax)
        if i == 0:
            ax.set_ylabel("Phase 2 / Phase 1 sites")
            ax.text(spec["xlim"][1], P2_N / P1_N * 1.08, f"{P2_N / P1_N:.0f}× participants", ha="right", va="bottom", fontsize=5.6, color="#555555")
        else:
            ax.tick_params(labelleft=False)
    x0, x1 = axes[0].get_position().x0, axes[-1].get_position().x1
    fig.text((x0 + x1) / 2, axes[0].get_position().y0 - 0.032, "Variant length (bp)", ha="center", va="top")
    return axes[0]


def render() -> tuple[Path, Path]:
    side = 3.80
    fig_p = new_figure(side, side)
    ax_p = fig_p.add_axes([0.18, 0.16, 0.76, 0.76])
    _panel_counts(ax_p, title=False)
    PANELDIR.mkdir(parents=True, exist_ok=True)
    fig_p.savefig(PANELDIR / "figS18_site_counts.pdf", dpi=600, bbox_inches=None)
    fig_p.savefig(PANELDIR / "figS18_site_counts.png", dpi=400, bbox_inches=None)
    plt.close(fig_p)

    fig = new_figure(9.3)
    gs = GridSpec(
        4,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.962,
        bottom=0.06,
        wspace=0.32,
        hspace=0.62,
        height_ratios=[0.42, 1.15, 1.12, 0.75],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_e = fig.add_subplot(gs[2, 1])

    _panel_scale(ax_a)
    label_panel(ax_a, "A", x=-0.01, y=1.12)
    _panel_counts(ax_b)
    label_panel(ax_b, "B", x=-0.12, y=1.06)
    _panel_svlen(ax_c)
    label_panel(ax_c, "C", x=-0.14, y=1.06)
    _panel_context(ax_d)
    label_panel(ax_d, "D", x=-0.12, y=1.06)
    _panel_extra(ax_e)
    label_panel(ax_e, "E", x=-0.12, y=1.06)
    ax_f = _panel_fold_by_length(fig, gs[3, :])
    label_panel(ax_f, "F", x=-0.12, y=1.12)

    stamp_mockup(fig, extra="Counts from Table 2")
    return save(fig, "figS18_p1p2_counts")
