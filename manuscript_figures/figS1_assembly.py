"""Figure S1: Sequencing coverage, diploid assemblies, and pedigrees."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyBboxPatch

from .style import (
    COL_IN,
    OKABE,
    P2_HIGH_N,
    P2_MID_N,
    P2_ONT_N,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)


def _panel_coverage(ax, rng: np.random.Generator) -> None:
    mid = rng.normal(16.3, 3.6, P2_MID_N)
    high = rng.normal(33.2, 4.3, P2_HIGH_N)
    ont = rng.normal(34.5, 5.1, P2_ONT_N)
    mid = mid[(mid > 6) & (mid < 32)]
    high = high[(high > 18) & (high < 50)]
    ont = ont[(ont > 15) & (ont < 55)]
    bins = np.linspace(5, 52, 38)
    ax.hist(mid, bins=bins, color=OKABE["blue"], alpha=0.75, label=f"PacBio mid-pass  n = {P2_MID_N:,}", zorder=2)
    ax.hist(high, bins=bins, color=OKABE["orange"], alpha=0.80, label=f"PacBio high-pass  n = {P2_HIGH_N:,}", zorder=3)
    ax.hist(ont, bins=bins, color=OKABE["green"], alpha=0.45, histtype="step", linewidth=1.3, label=f"ONT overlap  n = {P2_ONT_N:,}")
    ax.axvline(16.3, color=OKABE["blue"], ls=":", lw=0.7)
    ax.axvline(33.2, color=OKABE["orange"], ls=":", lw=0.7)
    ax.set_xlabel("Mean coverage (×)")
    ax.set_ylabel("Participants")
    ax.set_title("Coverage by sequencing stratum", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.6)


def _panel_ngx(ax) -> None:
    x = np.linspace(0.02, 0.98, 80)
    # Schematic NGx: high-pass remains long well past 50%; mid-pass drops after ~10 Mb.
    mid = 42 * np.exp(-3.4 * x) + 0.4
    high = 118 * np.exp(-1.15 * x) + 1.2
    ax.plot(x * 100, mid, color=OKABE["blue"], lw=1.5, label="Mid-pass  NG50 = 10.1 Mb")
    ax.plot(x * 100, high, color=OKABE["orange"], lw=1.5, label="High-pass  NG50 = 77.0 Mb")
    ax.axvline(50, color="#CCCCCC", ls="--", lw=0.6)
    ax.set_xlabel("Genome fraction  x  (%)")
    ax.set_ylabel("Contig length (Mb)")
    ax.set_title("Diploid assembly NGx (placeholder curves)", loc="left", pad=2)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 125)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.6)


def _family_glyph(ax, x0: float, y0: float, kind: str, color: str) -> None:
    """Tiny pedigree icon: trio, quartet, or three-generation."""
    r = 0.11

    def person(x, y, filled=True):
        ax.add_patch(
            Circle(
                (x, y),
                r,
                facecolor=color if filled else "white",
                edgecolor=color,
                linewidth=0.8,
                zorder=3,
                transform=ax.transData,
                clip_on=False,
            )
        )

    if kind == "trio":
        person(x0 - 0.22, y0 + 0.28)
        person(x0 + 0.22, y0 + 0.28)
        ax.plot([x0 - 0.22, x0 + 0.22], [y0 + 0.28, y0 + 0.28], color=color, lw=0.7)
        ax.plot([x0, x0], [y0 + 0.28, y0 - 0.05], color=color, lw=0.7)
        person(x0, y0 - 0.16)
    elif kind == "quartet":
        person(x0 - 0.22, y0 + 0.28)
        person(x0 + 0.22, y0 + 0.28)
        ax.plot([x0 - 0.22, x0 + 0.22], [y0 + 0.28, y0 + 0.28], color=color, lw=0.7)
        ax.plot([x0, x0], [y0 + 0.28, y0 + 0.02], color=color, lw=0.7)
        ax.plot([x0 - 0.18, x0 + 0.18], [y0 + 0.02, y0 + 0.02], color=color, lw=0.7)
        person(x0 - 0.18, y0 - 0.16)
        person(x0 + 0.18, y0 - 0.16)
    else:
        person(x0 - 0.28, y0 + 0.38)
        person(x0 + 0.04, y0 + 0.38)
        ax.plot([x0 - 0.28, x0 + 0.04], [y0 + 0.38, y0 + 0.38], color=color, lw=0.7)
        ax.plot([x0 - 0.12, x0 - 0.12], [y0 + 0.38, y0 + 0.10], color=color, lw=0.7)
        person(x0 - 0.28, y0 + 0.0)
        person(x0 + 0.04, y0 + 0.0)
        ax.plot([x0 - 0.28, x0 + 0.04], [y0 + 0.0, y0 + 0.0], color=color, lw=0.7)
        ax.plot([x0 - 0.12, x0 - 0.12], [y0 + 0.0, y0 - 0.28], color=color, lw=0.7)
        person(x0 - 0.12, y0 - 0.40)


def _panel_pedigrees(ax) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.2)
    ax.axis("off")
    rows = [
        (1.3, "trio", 25, 75, "Trios", OKABE["blue"]),
        (5.0, "quartet", 12, 48, "Quartets", OKABE["orange"]),
        (8.5, "multigen", 2, 11, "Multi-generational", OKABE["green"]),
    ]
    ax.set_title("Complete pedigrees in the discovery set", loc="left", pad=2)
    for x, kind, n_fam, n_ppl, lab, color in rows:
        _family_glyph(ax, x, 2.55, kind, color)
        ax.text(x, 1.35, lab, ha="center", va="center", fontsize=6.4, color="#222222")
        ax.text(x, 0.85, f"{n_fam} families", ha="center", va="center", fontsize=6.2, color=color)
        ax.text(x, 0.42, f"{n_ppl} people", ha="center", va="center", fontsize=5.8, color="#555555")
    ax.text(
        0.15,
        3.85,
        "39 complete pedigrees  ·  134 people  ·  32 wholly mid-pass, 2 high-pass, 5 mixed",
        fontsize=6.0,
        color="#444444",
        ha="left",
    )


def _panel_assembly_metrics(ax) -> None:
    labels = [
        "Assembled\nbases (Gbp)",
        "Contigs\n(thousands)",
        "Diploid NG50\n(Mb)",
        "Phased hap.\nNG50 (Mb)",
    ]
    mid = np.array([5.94, 7.823, 10.1, 2.0])
    high = np.array([6.05, 0.725, 77.0, 50.3])
    # Scale each pair so the larger value = 1 for a compact grouped view, annotate raw.
    scale = np.maximum(mid, high)
    x = np.arange(len(labels))
    w = 0.36
    ax.bar(x - w / 2, mid / scale, width=w, color=OKABE["blue"], label="Mid-pass  n = 11,128", zorder=2)
    ax.bar(x + w / 2, high / scale, width=w, color=OKABE["orange"], label="High-pass  n = 1,133", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Relative to stratum maximum")
    ax.set_title("Assembly continuity (Table 2)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.6)
    raw_mid = ["5.94", "7,823", "10.1", "2.0"]
    raw_high = ["6.05", "725", "77.0", "50.3"]
    for i, (a, b, sm) in enumerate(zip(raw_mid, raw_high, scale)):
        ax.text(i - w / 2, mid[i] / sm + 0.03, a, ha="center", va="bottom", fontsize=5.3, color=OKABE["blue"])
        ax.text(i + w / 2, high[i] / sm + 0.03, b, ha="center", va="bottom", fontsize=5.3, color=OKABE["orange"])
    ax.set_ylim(0, 1.22)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(4)
    ax_panel("figS1_coverage", 2.85, lambda ax: _panel_coverage(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS1_ngx", 2.85, _panel_ngx, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS1_pedigrees", 2.85, _panel_pedigrees, width=COL_IN * 1.12, left=0.08, right=0.96, top=0.90, bottom=0.08)
    ax_panel("figS1_assembly_metrics", 2.85, _panel_assembly_metrics, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    fig = new_figure(7.35)
    gs = GridSpec(
        2,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.95,
        bottom=0.07,
        wspace=0.32,
        hspace=0.42,
        height_ratios=[1.05, 1.0],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    _panel_coverage(ax_a, rng)
    _panel_ngx(ax_b)
    _panel_pedigrees(ax_c)
    _panel_assembly_metrics(ax_d)

    stamp_mockup(fig, extra="Coverage histograms are simulated around Table 1 means")
    return save(fig, "figS1_assembly")
