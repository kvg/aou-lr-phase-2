"""Figure S4: Pathogenic tandem-repeat loci resolved by TRGT."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle

from .style import (
    ANC_COLORS,
    COL_IN,
    OKABE,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)


def _panel_schematic(ax) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.6)
    ax.axis("off")
    ax.set_title("HiFi spanning of pathogenic tandem repeats", loc="left", pad=2)
    # Read as a long bar covering the repeat.
    ax.add_patch(Rectangle((0.4, 2.15), 9.1, 0.55, facecolor="#E8F1F8", edgecolor=OKABE["blue"], linewidth=0.8))
    ax.text(5.0, 2.42, "PacBio HiFi read  (~15–20 kb)", ha="center", va="center", fontsize=6.2, color=OKABE["blue"])
    parts = [
        (0.7, 1.6, "#D9D9D9", "flank"),
        (2.3, 3.4, OKABE["orange"], "repeat (motif × n)"),
        (5.7, 1.6, "#D9D9D9", "flank"),
    ]
    for x, w, color, lab in parts:
        ax.add_patch(Rectangle((x, 1.15), w, 0.55, facecolor=color, edgecolor="none", alpha=0.9))
        ax.text(x + w / 2, 1.42, lab, ha="center", va="center", fontsize=5.8, color="#222222")
    ax.annotate("", xy=(7.3, 2.15), xytext=(7.3, 1.70), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.8, mutation_scale=7))
    ax.text(
        0.4,
        0.45,
        "Primary association variable: total repeat-unit count   ·   Covariates: motif composition, longest pure segment",
        fontsize=5.8,
        color="#444444",
        ha="left",
    )
    ax.text(
        0.4,
        0.12,
        "Loci: FMR1, HTT, FXN, C9ORF72, RFC1, FAME  (full distributions in the companion repeat manuscript)",
        fontsize=5.6,
        color="#666666",
        ha="left",
    )


def _panel_lengths(ax, rng: np.random.Generator, locus: str, unit: str, bins, thresholds, ancs) -> None:
    xs = np.linspace(bins[0], bins[1], 240)
    specs = {
        "FMR1": [("AFR", 29, 4.0), ("EUR", 30, 3.5), ("EAS", 28, 3.2), ("AMR", 30, 4.2)],
        "HTT": [("AFR", 17, 3.2), ("EUR", 18, 2.8), ("EAS", 16, 2.4), ("AMR", 18, 3.0)],
    }[locus]
    for anc, mu, sd in specs:
        ys = np.exp(-0.5 * ((xs - mu) / sd) ** 2)
        # Long tail of intermediate / premutation alleles (placeholder).
        tail = 0.12 * np.exp(-0.5 * ((xs - (mu + 2.4 * sd)) / (sd * 1.6)) ** 2)
        ys = (ys + tail) / (ys + tail).max()
        ax.plot(xs, ys, color=ANC_COLORS[anc], lw=1.15, label=anc)
        ax.fill_between(xs, 0, ys, color=ANC_COLORS[anc], alpha=0.10, linewidth=0)
    colors_th = [OKABE["green"], OKABE["orange"], OKABE["vermillion"]]
    lifts = [1.02, 1.12, 1.02]
    x_nudge = [0.0, 0.0, 1.6]
    for i, (name, x) in enumerate(thresholds.items()):
        ax.axvline(x, color=colors_th[i], ls="--", lw=0.8)
        lift = lifts[i] if i < len(lifts) else 1.02
        dx = x_nudge[i] if i < len(x_nudge) else 0.0
        ax.text(x + dx, lift, name, ha="center", va="bottom", fontsize=5.2, color=colors_th[i])
    ax.set_xlim(bins[0], bins[1])
    ax.set_ylim(0, 1.18)
    ax.set_yticks([])
    ax.set_xlabel(f"{locus}  {unit} (placeholder)")
    ax.set_title(f"{locus} allele-length density", loc="left", pad=2)
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)
    ax.legend(loc="upper right", ncol=2, fontsize=5.5)


def _panel_concordance(ax, rng: np.random.Generator) -> None:
    n = 280
    trgt = rng.lognormal(np.log(28), 0.28, n)
    eh = trgt * rng.uniform(0.82, 1.05, n) + rng.normal(0, 1.8, n)
    eh = np.clip(eh, 5, None)
    # Short-read dropout on long alleles.
    long = trgt > 55
    eh[long] = np.clip(eh[long] * rng.uniform(0.45, 0.75, long.sum()), 20, 80)
    ax.scatter(trgt[~long], eh[~long], s=8, c=OKABE["blue"], alpha=0.45, linewidths=0, rasterized=True, label="Concordant range")
    ax.scatter(trgt[long], eh[long], s=10, c=OKABE["vermillion"], alpha=0.7, linewidths=0, rasterized=True, label="EH under-calls long alleles")
    ax.plot([5, 90], [5, 90], ls="--", color="#888888", lw=0.7)
    ax.set_xlim(5, 95)
    ax.set_ylim(5, 95)
    ax.set_xlabel("TRGT length (HiFi, overlapping samples)")
    ax.set_ylabel("ExpansionHunter length (srWGS)")
    ax.set_title("TRGT vs ExpansionHunter (placeholder)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left", fontsize=5.5)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(14)
    ax_panel("figS4_schematic", 1.85, _panel_schematic, left=0.03, right=0.99, top=0.86, bottom=0.06)

    def _fmr1(ax):
        _panel_lengths(ax, rng, "FMR1", "CGG copies", (10, 90), {"normal": 45, "premut.": 55}, None)
        ax.set_xlim(10, 90)
        ax.text(0.98, 0.90, "full >200 CGG off-scale", transform=ax.transAxes, ha="right", fontsize=5.4, color=OKABE["vermillion"])

    ax_panel("figS4_fmr1", 2.85, _fmr1, width=COL_IN * 1.12, left=0.12, right=0.96, top=0.88, bottom=0.16)
    ax_panel(
        "figS4_htt",
        2.85,
        lambda ax: _panel_lengths(ax, rng, "HTT", "CAG copies", (8, 50), {"normal": 27, "reduced": 36, "full": 40}, None),
        width=COL_IN * 1.12,
        left=0.12,
        right=0.96,
        top=0.88,
        bottom=0.16,
    )
    ax_panel("figS4_concordance", 2.55, lambda ax: _panel_concordance(ax, rng), left=0.10, right=0.98, top=0.88, bottom=0.16)
    fig = new_figure(7.50)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.95,
        bottom=0.07,
        wspace=0.30,
        hspace=0.50,
        height_ratios=[0.78, 1.12, 1.12],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])

    _panel_schematic(ax_a)
    _panel_lengths(
        ax_b,
        rng,
        "FMR1",
        "CGG copies",
        (10, 90),
        {"normal": 45, "premut.": 55},
        None,
    )
    # full-mutation threshold is off-axis; annotate instead of drawing x=200.
    ax_b.set_xlim(10, 90)
    ax_b.text(0.98, 0.90, "full >200 CGG off-scale", transform=ax_b.transAxes, ha="right", fontsize=5.4, color=OKABE["vermillion"])
    _panel_lengths(
        ax_c,
        rng,
        "HTT",
        "CAG copies",
        (8, 50),
        {"normal": 27, "reduced": 36, "full": 40},
        None,
    )
    _panel_concordance(ax_d, rng)

    stamp_mockup(fig, extra="Allele lengths and EH comparison are layout placeholders")
    return save(fig, "figS4_repeats")
