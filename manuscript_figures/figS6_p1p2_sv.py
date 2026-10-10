"""Figure S6: Phase 1 vs Phase 2 SV catalog comparison (placeholder overlap)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from .style import (
    ANC_COLORS,
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
    subplot_panel,
)

# Placeholder site counts (≥50 bp). Phase 1 catalog ~650k from the imputation section.
P1_SV = 650_000
P2_SV = SV_N["SV_50"]
SHARED = 580_000
P1_ONLY = P1_SV - SHARED
P2_ONLY = P2_SV - SHARED

AF_CMAP = LinearSegmentedColormap.from_list("af", ["#F4F8FB", OKABE["blue"]])


def _panel_accounting(ax) -> None:
    """Two-row stacked bars: each catalog split into shared / unique."""
    ax.set_xlim(0, P2_SV * 1.16)
    ax.set_ylim(-0.55, 2.45)
    ax.axis("off")
    ax.set_title("SV ≥50 bp site accounting (placeholder Truvari overlap)", loc="left", pad=4)
    rows = [
        (1.20, "Phase 1", P1_SV, SHARED, P1_ONLY, OKABE["blue"], OKABE["vermillion"]),
        (0.20, "Phase 2", P2_SV, SHARED, P2_ONLY, OKABE["blue"], OKABE["orange"]),
    ]
    for y, name, total, shared, unique, c_sh, c_un in rows:
        ax.add_patch(Rectangle((0, y), shared, 0.72, facecolor=c_sh, edgecolor="none", zorder=2))
        ax.add_patch(Rectangle((shared, y), unique, 0.72, facecolor=c_un, edgecolor="none", zorder=2))
        ax.text(-0.012 * P2_SV, y + 0.36, name, ha="right", va="center", fontsize=7.0, color="#222222")
        ax.text(shared / 2, y + 0.36, f"shared  {format_millions(shared)}", ha="center", va="center", fontsize=6.0, color="white")
        if unique > 0.10 * P2_SV:
            ax.text(
                shared + unique / 2,
                y + 0.36,
                f"unique  {format_millions(unique)}",
                ha="center",
                va="center",
                fontsize=6.0,
                color="#222222",
            )
        ax.text(total + 0.015 * P2_SV, y + 0.50, format_millions(total), ha="left", va="center", fontsize=6.2, color="#333333")
        if unique <= 0.10 * P2_SV:
            ax.text(
                total + 0.015 * P2_SV,
                y + 0.18,
                f"{format_millions(unique)} unique",
                ha="left",
                va="center",
                fontsize=5.6,
                color=OKABE["vermillion"],
            )
    ax.text(
        0,
        2.12,
        f"Almost all Phase 1 sites recovered  ·  {P2_ONLY / P2_SV:.0%} of Phase 2 is new  ·  n = {P1_N:,} vs {P2_N:,}",
        fontsize=6.0,
        color="#444444",
    )


def _panel_upset(ax):
    """Three-way overlap: HPRC, Phase 2, Phase 1 (row order)."""
    # Keys are (in_HPRC, in_P2, in_P1). Counts are layout placeholders.
    combos = {
        (0, 1, 0): 1_150_000,  # P2 private
        (1, 1, 0): 770_000,  # P2 ∩ HPRC, not P1
        (1, 1, 1): 400_000,  # all three
        (0, 1, 1): 180_000,  # P1 ∩ P2, not HPRC
        (1, 0, 0): 210_000,  # HPRC only
        (0, 0, 1): 50_000,  # P1 only
        (1, 0, 1): 20_000,  # P1 ∩ HPRC, lost in P2
    }
    keys = sorted(combos, key=combos.get, reverse=True)
    counts = np.array([combos[k] for k in keys], dtype=float)
    n = len(keys)
    fig = ax.figure
    spec = ax.get_subplotspec()
    ax.remove()
    inner = spec.subgridspec(2, 2, width_ratios=[0.30, 0.70], height_ratios=[0.62, 0.38], wspace=0.08, hspace=0.10)
    ax_bar = fig.add_subplot(inner[0, 1])
    ax_mat = fig.add_subplot(inner[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(inner[1, 0])

    ax_bar.bar(np.arange(n), counts / 1e6, color="#3D3D3D", width=0.78, zorder=2)
    ax_bar.set_ylabel("Sites (millions)")
    ax_bar.tick_params(labelbottom=False)
    despine(ax_bar)
    ax_bar.set_title("Phase 1 × Phase 2 × HPRC (placeholder)", loc="left", pad=2)
    for i, c in enumerate(counts):
        ax_bar.text(i, c / 1e6 + 0.04, format_millions(c), ha="center", va="bottom", fontsize=5.2)

    labels = ["HPRC", "Phase 2", "Phase 1"]
    mat = np.array(keys, dtype=float).T
    ax_mat.set_xlim(-0.5, n - 0.5)
    ax_mat.set_ylim(-0.5, 2.5)
    ax_mat.set_yticks([])
    ax_mat.set_xticks([])
    for r in range(3):
        ax_mat.plot([-0.4, n - 0.6], [r, r], color="#EFEFEF", lw=6, zorder=1)
        for c in range(n):
            on = mat[r, c] > 0
            ax_mat.scatter(c, r, s=26, c="#222222" if on else "#E2E2E2", zorder=3, linewidths=0)
    for c, key in enumerate(keys):
        ys = [i for i, v in enumerate(key) if v]
        if len(ys) >= 2:
            ax_mat.plot([c, c], [min(ys), max(ys)], color="#222222", lw=0.85, zorder=2)
    despine(ax_mat, left=True, bottom=False)
    ax_mat.tick_params(bottom=False)

    set_n = [sum(combos[k] for k in keys if k[i]) for i in range(3)]
    ax_set.barh([0, 1, 2], np.array(set_n) / 1e6, color="#4C4C4C", height=0.55)
    ax_set.set_yticks([0, 1, 2])
    ax_set.set_yticklabels(labels, fontsize=6.0)
    ax_set.yaxis.tick_right()
    ax_set.tick_params(axis="y", length=0, pad=3)
    ax_set.set_ylim(-0.5, 2.5)
    ax_set.invert_xaxis()
    despine(ax_set, left=False)
    ax_set.spines["right"].set_visible(False)
    ax_set.set_xticks([])
    return ax_bar


def _panel_private_source(ax) -> None:
    """Phase 2–only sites: HPRC-known vs ancestry-private."""
    labels = [
        "In HPRC,\nabsent from P1",
        "AFR-private\n(more n, 15–30×)",
        "AMR-private",
        "EAS-private",
        "EUR-private",
        "MID/SAS-\nprivate",
    ]
    vals = np.array([0.40, 0.28, 0.14, 0.07, 0.05, 0.06]) * P2_ONLY
    colors = [
        "#B8B8B8",
        ANC_COLORS["AFR"],
        ANC_COLORS["AMR"],
        ANC_COLORS["EAS"],
        ANC_COLORS["EUR"],
        ANC_COLORS["SAS"],
    ]
    y = np.arange(len(labels))
    ax.barh(y, vals / 1e6, color=colors, height=0.72, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("Phase 2–only sites (millions)")
    ax.set_title("Where Phase 2–only SVs come from", loc="left", pad=2)
    despine(ax)
    for yi, v in zip(y, vals):
        ax.text(
            v / 1e6 + 0.02,
            yi,
            f"{v / P2_ONLY:.0%}  ·  {format_millions(v)}",
            va="center",
            fontsize=5.5,
            color="#333333",
        )
    ax.set_xlim(0, vals.max() / 1e6 * 1.38)


def _panel_af(ax, rng: np.random.Generator) -> None:
    n = 8000
    af2 = rng.beta(0.45, 8.0, n)
    af1 = np.clip(af2 * rng.uniform(0.75, 1.25, n) + rng.normal(0, 0.004, n), 1e-4, 0.95)
    # AFR-only Phase 1 overestimates AF for variants that are rare in the pan-ancestry catalog.
    rare = af2 < 0.02
    af1[rare] = np.clip(af1[rare] * rng.uniform(1.1, 2.4, rare.sum()), 1e-4, 0.2)
    hb = ax.hexbin(
        af1,
        af2,
        gridsize=38,
        xscale="log",
        yscale="log",
        cmap=AF_CMAP,
        mincnt=2,
        linewidths=0,
    )
    ax.plot([1e-4, 1], [1e-4, 1], ls="--", color="#888888", lw=0.7)
    ax.set_xlim(8e-4, 1)
    ax.set_ylim(8e-4, 1)
    ax.set_xlabel("Phase 1 allele frequency")
    ax.set_ylabel("Phase 2 allele frequency")
    ax.set_title("Shared sites, AF vs AF", loc="left", pad=2)
    despine(ax)
    ax.text(0.05, 0.92, "shared n ≈ 0.58M (placeholder)", transform=ax.transAxes, fontsize=5.6, color="#444444")
    cb = ax.figure.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=5.2)
    cb.set_label("sites", fontsize=5.5)


def _panel_discovery(ax, rng: np.random.Generator) -> None:
    """Unique SVs vs n: AFR-only saturation vs full pan-ancestry cohort."""
    n_afr = np.array([50, 100, 200, 400, 800, 1027, 1500, 2000, 2816])
    n_all = np.array([200, 500, 1027, 2000, 4000, 7000, 10000, 12261])

    def saturating(n, k, xmax):
        return k * (1.0 - np.exp(-n / xmax)) + rng.normal(0, 0.008, n.size) * k

    afr = saturating(n_afr, 1.05e6, 2200)
    full = saturating(n_all, 2.55e6, 5500)
    afr[-1] = 1.12e6
    full[-1] = P2_SV
    ax.plot(n_afr, afr / 1e6, color=ANC_COLORS["AFR"], lw=1.4, marker="o", ms=3.4, label="AFR-only samples")
    ax.plot(n_all, full / 1e6, color=OKABE["blue"], lw=1.4, marker="o", ms=3.4, label="Pan-ancestry cohort")
    ax.axvline(P1_N, color="#AAAAAA", ls=":", lw=0.7)
    ax.axvline(P2_N, color="#AAAAAA", ls=":", lw=0.7)
    ax.text(P1_N + 80, 0.22, "Phase 1 n", fontsize=5.4, color="#666666", rotation=90, va="bottom")
    ax.text(P2_N - 380, 0.22, "Phase 2 n", fontsize=5.4, color="#666666", rotation=90, va="bottom")
    ax.set_xlabel("Discovery samples")
    ax.set_ylabel("Unique SVs ≥50 bp (millions)")
    ax.set_title("Discovery curves (placeholder)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="lower right", fontsize=5.5)
    ax.set_xlim(0, 13000)
    ax.set_ylim(0, 2.85)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(6)
    ax_panel("figS6_accounting", 1.85, _panel_accounting, left=0.06, right=0.98, top=0.86, bottom=0.10)
    subplot_panel("figS6_upset", 3.25, _panel_upset, width=COL_IN * 1.25, left=0.14, right=0.96, top=0.90, bottom=0.08)
    ax_panel("figS6_private_source", 2.85, _panel_private_source, width=COL_IN * 1.12, left=0.22, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS6_af", 2.85, lambda ax: _panel_af(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS6_discovery", 2.85, lambda ax: _panel_discovery(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    fig = new_figure(7.50)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.10,
        right=0.98,
        top=0.94,
        bottom=0.06,
        wspace=0.32,
        hspace=0.50,
        height_ratios=[0.72, 1.12, 1.12],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b_host = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_e = fig.add_subplot(gs[2, 1])

    _panel_accounting(ax_a)
    _panel_upset(ax_b_host)
    _panel_private_source(ax_c)
    _panel_af(ax_d, rng)
    _panel_discovery(ax_e, rng)

    stamp_mockup(fig, extra="Overlap counts are layout placeholders pending Truvari of P1 vs P2 vs HPRC")
    return save(fig, "figS6_p1p2_sv")
