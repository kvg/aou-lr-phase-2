"""Figure S10: AoU novelty relative to HPRC ∪ HGSVC3, by class and context."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    P2_ANC_N,
    P2_N,
    SV_N,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    stamp_mockup,
)

C_KNOWN = "#C8C8C8"
C_NOVEL = OKABE["blue"]

# Placeholder Truvari splits of ≥50 bp INS/DEL vs HPRC ∪ HGSVC3 dipcall catalogs.
# Columns: title, n_aou, n_panel, n_known (AoU ∩ panel), tau_known, tau_novel, recall_inf
FACETS = (
    ("INS, unique sequence", 0.22 * SV_N["INS_50"], 210_000, 0.42, 380, 4200, 0.91),
    ("INS, repetitive", 0.78 * SV_N["INS_50"], 780_000, 0.35, 900, 6500, 0.72),
    ("DEL, unique sequence", 0.38 * SV_N["DEL_50"], 165_000, 0.58, 320, 3800, 0.93),
    ("DEL, repetitive", 0.62 * SV_N["DEL_50"], 280_000, 0.45, 750, 5200, 0.76),
)


def _grid_x(n: int) -> np.ndarray:
    early = np.linspace(1, 250, 36)
    late = np.geomspace(250, n, 70)
    return np.unique(np.concatenate([early, late])).astype(int)


def _rarefy(x: np.ndarray, n_aou: float, frac_known: float, tau_k: float, tau_n: float):
    known_inf = n_aou * frac_known
    novel_inf = n_aou * (1.0 - frac_known)
    known = known_inf * (1.0 - np.exp(-x / tau_k))
    novel = novel_inf * (1.0 - np.exp(-x / tau_n))
    total = known + novel
    iqr = 0.07 * total / np.sqrt(np.clip(x, 8, None) / 40.0)
    return known, novel, total, iqr


def _panel_rarefy(ax, spec: tuple, x: np.ndarray, ylabel: bool, legend: bool) -> None:
    title, n_aou, n_panel, frac_known, tau_k, tau_n, _rec = spec
    known, novel, total, iqr = _rarefy(x, n_aou, frac_known, tau_k, tau_n)
    ax.fill_between(x, 0, known, color=C_KNOWN, lw=0, label="Also in HPRC + HGSVC3")
    ax.fill_between(x, known, total, color=C_NOVEL, lw=0, alpha=0.88, label="AoU-only (Truvari miss)")
    ax.fill_between(x, np.clip(total - iqr, 0, None), total + iqr, color="#111111", alpha=0.10, lw=0)
    ax.axhline(n_panel, color="#666666", ls=":", lw=0.9)
    ax.text(x[-1] * 0.98, n_panel, "HPRC+HGSVC3", ha="right", va="bottom", fontsize=5.3, color="#555555")
    ax.set_xlim(1, P2_N)
    ax.set_title(title, loc="left", pad=2)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: format_millions(v) if v >= 1000 else f"{v:.0f}"))
    if ylabel:
        ax.set_ylabel("Unique SVs ≥50 bp")
    despine(ax)
    ax.text(0.04, 0.90, format_millions(n_aou), transform=ax.transAxes, fontsize=6.0, color="#333333")
    if legend:
        ax.legend(loc="upper left", fontsize=5.3, bbox_to_anchor=(0.0, 0.86))


def _panel_recall(ax, x: np.ndarray) -> None:
    styles = [
        (OKABE["blue"], "-", "INS unique"),
        (OKABE["orange"], "-", "INS repetitive"),
        (OKABE["blue"], "--", "DEL unique"),
        (OKABE["orange"], "--", "DEL repetitive"),
    ]
    for spec, (color, ls, lab) in zip(FACETS, styles):
        _, _, n_panel, _, tau_k, _, rec_inf = spec
        # Recall of the *panel* catalog, not of AoU.
        y = rec_inf * (1.0 - np.exp(-x / (tau_k * 1.15)))
        ax.plot(x, y, color=color, ls=ls, lw=1.35, label=lab)
    ax.set_xlim(1, P2_N)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Fraction of panel sites recovered")
    ax.set_xlabel("AoU samples (random order)")
    ax.set_title("Recall of HPRC + HGSVC3", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="lower right", fontsize=5.3)


def _panel_per_genome(ax, rng: np.random.Generator) -> None:
    # Novel SVs per genome vs the combined assembly panel (placeholder).
    # Unique vs repetitive, by ancestry. AFR highest leftover.
    means_u = {"AFR": 1450, "AMR": 980, "SAS": 820, "MID": 760, "EAS": 710, "OTH": 900, "EUR": 620}
    means_r = {"AFR": 3100, "AMR": 2100, "SAS": 1750, "MID": 1600, "EAS": 1480, "OTH": 1900, "EUR": 1320}
    order = [a for a in ANC_ORDER if a in P2_ANC_N]
    x = np.arange(len(order))
    w = 0.38
    yu = np.array([means_u[a] for a in order], dtype=float)
    yr = np.array([means_r[a] for a in order], dtype=float)
    yu *= rng.uniform(0.97, 1.03, yu.size)
    yr *= rng.uniform(0.97, 1.03, yr.size)
    ax.bar(x - w / 2, yu, width=w, color=OKABE["blue"], label="Unique sequence", zorder=2)
    ax.bar(x + w / 2, yr, width=w, color=OKABE["orange"], label="Repetitive", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("AoU-only SVs / genome")
    ax.set_title("Per-genome leftover after the assembly panels", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.3)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(23)
    x = _grid_x(P2_N)
    slugs = ["ins_unique", "ins_repetitive", "del_unique", "del_repetitive"]
    ymax_ins = max(FACETS[0][1], FACETS[1][1]) * 1.08
    ymax_del = max(FACETS[2][1], FACETS[3][1]) * 1.08
    for i, (spec, slug) in enumerate(zip(FACETS, slugs)):
        def _facet(ax, spec=spec, i=i):
            _panel_rarefy(ax, spec, x, ylabel=True, legend=(i == 0))
            ax.set_ylim(0, ymax_ins if i < 2 else ymax_del)
            ax.set_xlabel("AoU samples (mean of random orders)")

        ax_panel(f"figS10_{slug}", 2.65, _facet, left=0.12, right=0.98, top=0.86, bottom=0.16)
    ax_panel("figS10_recall", 2.65, lambda ax: _panel_recall(ax, x), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS10_per_genome", 2.65, lambda ax: _panel_per_genome(ax, rng), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)

    fig = new_figure(7.85)
    gs = GridSpec(
        2,
        1,
        figure=fig,
        left=0.09,
        right=0.98,
        top=0.93,
        bottom=0.07,
        hspace=0.28,
        height_ratios=[1.55, 1.0],
    )
    top = gs[0].subgridspec(2, 2, wspace=0.26, hspace=0.42)
    bot = gs[1].subgridspec(1, 2, wspace=0.30)
    axes = [fig.add_subplot(top[i, j]) for i in range(2) for j in range(2)]
    ax_e = fig.add_subplot(bot[0, 0])
    ax_f = fig.add_subplot(bot[0, 1])

    for i, (ax, spec) in enumerate(zip(axes, FACETS)):
        _panel_rarefy(ax, spec, x, ylabel=(i % 2 == 0), legend=(i == 0))
        if i < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("AoU samples (mean of random orders)")

    ymax_ins = max(FACETS[0][1], FACETS[1][1]) * 1.08
    ymax_del = max(FACETS[2][1], FACETS[3][1]) * 1.08
    for ax in axes[:2]:
        ax.set_ylim(0, ymax_ins)
    for ax in axes[2:]:
        ax.set_ylim(0, ymax_del)

    _panel_recall(ax_e, x)
    _panel_per_genome(ax_f, rng)

    stamp_mockup(
        fig,
        extra="Bands are IQR across random sample orders. Panel overlap is a Truvari placeholder",
    )
    return save(fig, "figS10_panel_novelty")
