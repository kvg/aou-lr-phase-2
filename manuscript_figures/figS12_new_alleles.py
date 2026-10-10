"""Figure S12: New alleles at Phase 1 sites, not only new sites."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from .style import (
    COL_IN,
    OKABE,
    P1_N,
    P2_N,
    SV_N,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    stamp_mockup,
)

P1_SITES = 650_000
P2_SITES = SV_N["SV_50"]
SHARED = 580_000
P2_ONLY = P2_SITES - SHARED

# Distinct ALTs per matched site (placeholder). REF is not counted.
N_ALT_P1 = 1.18 * SHARED  # mostly biallelic
N_ALT_NEW_AT_OLD = 420_000  # extra ALTs Phase 2 adds at matched sites
N_ALT_ON_NEW = 1.08 * P2_ONLY  # new sites, still mostly one ALT
N_ALT_P2 = N_ALT_P1 + N_ALT_NEW_AT_OLD + N_ALT_ON_NEW

C_OLD = OKABE["blue"]
C_EXTRA = OKABE["green"]
C_NEW_SITE = OKABE["orange"]
HEAT = LinearSegmentedColormap.from_list("n", ["#F7F7F7", OKABE["blue"], "#001F3F"])


def _panel_accounting(ax) -> None:
    """Sites hide the extra-allele slice; allele accounting does not."""
    ax.set_xlim(0, N_ALT_P2 * 1.18)
    ax.set_ylim(-0.55, 2.55)
    ax.axis("off")
    ax.set_title("Phase 2 added sites and alleles at sites Phase 1 already had", loc="left", pad=4)

    # Sites row: two parts (shared / P2-only), scaled onto the allele x-axis for comparison.
    y_s, y_a = 1.35, 0.22
    ax.add_patch(Rectangle((0, y_s), SHARED, 0.72, facecolor=C_OLD, edgecolor="none", zorder=2))
    ax.add_patch(Rectangle((SHARED, y_s), P2_ONLY, 0.72, facecolor=C_NEW_SITE, edgecolor="none", zorder=2))
    ax.text(-0.01 * N_ALT_P2, y_s + 0.36, "Sites", ha="right", va="center", fontsize=7.0)
    ax.text(SHARED / 2, y_s + 0.36, f"matched to P1  {format_millions(SHARED)}", ha="center", va="center", fontsize=5.8, color="white")
    ax.text(SHARED + P2_ONLY / 2, y_s + 0.36, f"new sites  {format_millions(P2_ONLY)}", ha="center", va="center", fontsize=5.8, color="#222222")
    ax.text(P2_SITES + 0.012 * N_ALT_P2, y_s + 0.36, format_millions(P2_SITES), ha="left", va="center", fontsize=6.2)

    ax.add_patch(Rectangle((0, y_a), N_ALT_P1, 0.72, facecolor=C_OLD, edgecolor="none", zorder=2))
    ax.add_patch(Rectangle((N_ALT_P1, y_a), N_ALT_NEW_AT_OLD, 0.72, facecolor=C_EXTRA, edgecolor="none", zorder=2))
    ax.add_patch(Rectangle((N_ALT_P1 + N_ALT_NEW_AT_OLD, y_a), N_ALT_ON_NEW, 0.72, facecolor=C_NEW_SITE, edgecolor="none", zorder=2))
    ax.text(-0.01 * N_ALT_P2, y_a + 0.36, "Alleles", ha="right", va="center", fontsize=7.0)
    ax.text(N_ALT_P1 / 2, y_a + 0.36, "P1 ALTs\nrecovered", ha="center", va="center", fontsize=5.5, color="white")
    ax.text(
        N_ALT_P1 + N_ALT_NEW_AT_OLD / 2,
        y_a + 0.36,
        f"new ALTs at\nP1 sites  {format_millions(N_ALT_NEW_AT_OLD)}",
        ha="center",
        va="center",
        fontsize=5.5,
        color="white",
        fontweight="bold",
    )
    ax.text(
        N_ALT_P1 + N_ALT_NEW_AT_OLD + N_ALT_ON_NEW / 2,
        y_a + 0.36,
        f"ALTs on new sites  {format_millions(N_ALT_ON_NEW)}",
        ha="center",
        va="center",
        fontsize=5.5,
        color="#222222",
    )
    ax.text(N_ALT_P2 + 0.012 * N_ALT_P2, y_a + 0.36, format_millions(N_ALT_P2), ha="left", va="center", fontsize=6.2)
    ax.text(
        0,
        2.22,
        f"Green is invisible if you only count Truvari sites  ·  n = {P1_N:,} vs {P2_N:,}",
        fontsize=6.0,
        color="#444444",
    )


def _simulate_n_alt(rng: np.random.Generator, n: int, extra_mean: float) -> tuple[np.ndarray, np.ndarray]:
    n1 = rng.choice([1, 2, 3, 4, 5], size=n, p=[0.82, 0.11, 0.04, 0.02, 0.01])
    extra = rng.poisson(extra_mean, n)
    extra = np.where(n1 == 1, np.minimum(extra, rng.choice([0, 0, 0, 1, 2], n)), extra)
    n2 = np.clip(n1 + extra, 1, 24)
    return n1, n2


def _panel_cardinality(ax, rng: np.random.Generator) -> None:
    n_u, n_r = 8000, 8000
    u1, u2 = _simulate_n_alt(rng, n_u, extra_mean=0.18)
    r1, r2 = _simulate_n_alt(rng, n_r, extra_mean=1.35)
    # Weight repetitive more in the heatmap by stacking.
    n1 = np.concatenate([u1, r1])
    n2 = np.concatenate([u2, r2])
    maxn = 12
    H = np.zeros((maxn, maxn))
    for a, b in zip(n1, n2):
        if a <= maxn and b <= maxn:
            H[b - 1, a - 1] += 1
    # Rescale to catalog sizes (unique vs repetitive mix ~ 1:3 among shared).
    H *= SHARED / H.sum()
    im = ax.imshow(np.log10(H + 1), origin="lower", cmap=HEAT, aspect="equal", extent=(0.5, maxn + 0.5, 0.5, maxn + 0.5))
    ax.plot([0.5, maxn + 0.5], [0.5, maxn + 0.5], ls="--", color="#888888", lw=0.7)
    ax.set_xlabel("Distinct ALTs at site in Phase 1")
    ax.set_ylabel("Distinct ALTs at site in Phase 2")
    ax.set_title("Matched sites: allele cardinality", loc="left", pad=2)
    ax.set_xticks(range(1, maxn + 1, 2))
    ax.set_yticks(range(1, maxn + 1, 2))
    despine(ax)
    ax.text(0.06, 0.90, "above diagonal: new alleles at an old site", transform=ax.transAxes, fontsize=5.5, color="#444444")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("log10 matched sites", fontsize=5.5)
    cb.ax.tick_params(labelsize=5.2)


def _panel_where(ax) -> None:
    labels = ["INS unique", "INS repetitive", "DEL unique", "DEL repetitive"]
    extra = np.array([0.08, 0.62, 0.06, 0.24]) * N_ALT_NEW_AT_OLD
    colors = [OKABE["blue"], OKABE["orange"], OKABE["sky"], OKABE["vermillion"]]
    y = np.arange(len(labels))
    ax.barh(y, extra / 1e3, color=colors, height=0.68, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("New ALTs at Phase 1 sites (thousands)")
    ax.set_title("Extra alleles are mostly repetitive insertions", loc="left", pad=2)
    despine(ax)
    for yi, v in zip(y, extra):
        ax.text(v / 1e3 + 8, yi, f"{v / N_ALT_NEW_AT_OLD:.0%}  ·  {format_millions(v)}", va="center", fontsize=5.5, color="#333333")
    ax.set_xlim(0, extra.max() / 1e3 * 1.42)


def _panel_locus(ax, rng: np.random.Generator) -> None:
    """Placeholder VNTR: Phase 2 fills in lengths Phase 1 never saw."""
    p1 = rng.normal(18.0, 2.4, 1027)
    p1 = p1[(p1 > 10) & (p1 < 26)]
    p2_old = rng.normal(18.0, 2.6, 4000)
    p2_new = np.concatenate(
        [
            rng.normal(12.0, 1.1, 900),
            rng.normal(27.5, 1.6, 700),
            rng.normal(34.0, 2.0, 220),
        ]
    )
    p2 = np.concatenate([p2_old, p2_new])
    p2 = p2[(p2 > 6) & (p2 < 42)]
    bins = np.arange(7, 43, 1)
    ax.hist(p1, bins=bins, density=True, histtype="step", color=OKABE["blue"], lw=1.4, label="Phase 1 alleles")
    ax.hist(p2, bins=bins, density=True, histtype="stepfilled", color=OKABE["green"], alpha=0.35, label="Phase 2 alleles")
    ax.hist(p2, bins=bins, density=True, histtype="step", color=OKABE["green"], lw=1.15)
    ax.set_xlabel("Repeat copies at a placeholder VNTR (matched site)")
    ax.set_ylabel("Density")
    ax.set_title("Same site, longer allelic series", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.4)
    ax.annotate(
        "new lengths",
        xy=(34, 0.02),
        xytext=(28, 0.11),
        fontsize=5.6,
        color=OKABE["green"],
        arrowprops=dict(arrowstyle="-|>", color=OKABE["green"], lw=0.6, mutation_scale=6),
    )


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(12)
    ax_panel("figS12_accounting", 2.15, _panel_accounting, left=0.08, right=0.98, top=0.86, bottom=0.12)
    ax_panel("figS12_cardinality", 2.85, lambda ax: _panel_cardinality(ax, rng), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS12_where", 2.35, _panel_where, width=COL_IN * 1.15, left=0.18, right=0.96, top=0.86, bottom=0.16)
    ax_panel("figS12_locus", 2.35, lambda ax: _panel_locus(ax, rng), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.86, bottom=0.16)
    fig = new_figure(7.55)
    gs = GridSpec(
        2,
        2,
        figure=fig,
        left=0.10,
        right=0.98,
        top=0.93,
        bottom=0.08,
        wspace=0.32,
        hspace=0.38,
        height_ratios=[0.95, 1.15],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    # split bottom-right into where + locus
    right = gs[1, 1].subgridspec(2, 1, hspace=0.55)
    ax_c = fig.add_subplot(right[0])
    ax_d = fig.add_subplot(right[1])

    _panel_accounting(ax_a)
    _panel_cardinality(ax_b, rng)
    _panel_where(ax_c)
    _panel_locus(ax_d, rng)

    stamp_mockup(fig, extra="Allele = distinct ALT at a Truvari-matched site. Counts are layout placeholders")
    return save(fig, "figS12_new_alleles")
