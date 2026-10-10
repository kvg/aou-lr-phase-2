"""Figure S8: Ebert-style cumulative SV discovery (placeholder)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from .style import (
    ANC_COLORS,
    COL_IN,
    P1_N,
    P2_ANC_N,
    P2_N,
    SV_N,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    save_panel,
    stamp_mockup,
)

# Ebert / matplotlib-tab10 frequency colors (same as scripts/plot_discovery.py).
FREQ_ORDER = ("shared", "major", "polymorphic", "singleton")
FREQ_COLORS = {
    "shared": "#d62728",
    "major": "#9467bd",
    "polymorphic": "#1f77b4",
    "singleton": "#17becf",
}
FREQ_LABELS = {
    "shared": "Shared (all samples; absent from GRCh38)",
    "major": "Major (AF ≥ 50%)",
    "polymorphic": "Polymorphic (≥2 samples, AF < 50%)",
    "singleton": "Singleton",
}

# Locked discovery order: non-African first, then African (sv_site_utils.ancestry_sort_key).
# MID is unnamed in that key so it follows OTH and precedes AFR.
P2_ORDER = ("EAS", "AMR", "EUR", "SAS", "OTH", "MID", "AFR")


def _p2_blocks() -> list[tuple[str, int, int]]:
    blocks = []
    start = 0
    for anc in P2_ORDER:
        n = P2_ANC_N[anc]
        blocks.append((anc, start, start + n))
        start += n
    return blocks


def _incremental(
    n: int,
    per_genome: float,
    target: float,
    afr_start: int | None,
    afr_mult: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """New unique SVs contributed by sample k (0-based)."""
    new = np.empty(n, dtype=float)
    new[0] = per_genome
    idx = np.arange(1, n)
    raw = (idx + 1) ** -0.55
    if afr_start is not None:
        raw[idx >= afr_start] *= afr_mult
    raw *= rng.uniform(0.93, 1.07, raw.size)
    new[1:] = raw * (target - per_genome) / raw.sum()
    return new


def _strata(
    new: np.ndarray,
    shared_inf: float,
    major_plat: float,
    major_tau: float,
    sing_tau: float,
    sing_max: float,
    afr_start: int | None,
) -> dict[str, np.ndarray]:
    """Prefix-recomputed Ebert frequency classes (not the final-cohort labels)."""
    n = new.size
    k = np.arange(1, n + 1, dtype=float)
    total = np.cumsum(new)
    shared = shared_inf + (total[0] - shared_inf) * np.power(k, -0.82)
    if afr_start is not None:
        extra = np.ones(n)
        extra[afr_start:] = 1.0 / (1.0 + 0.0018 * (k[afr_start:] - afr_start))
        shared = shared * extra
        shared = np.minimum.accumulate(shared)
    shared = np.clip(shared, min(shared_inf, total[-1]), total)
    shared[0] = total[0]

    major = major_plat * (1.0 - np.exp(-(k - 1.0) / major_tau))
    major[0] = 0.0
    major = np.minimum(major, np.maximum(0.40 * total - shared, 0.0))
    overflow = shared + major - total
    major = np.where(overflow > 0, major - overflow, major)
    major = np.clip(major, 0, None)

    remain = np.maximum(total - shared - major, 0.0)
    frac = np.zeros(n)
    frac[1:] = sing_max * (1.0 - np.exp(-(k[1:] - 1.0) / sing_tau))
    singleton = remain * frac
    polymorphic = remain - singleton
    return {
        "shared": shared,
        "major": major,
        "polymorphic": polymorphic,
        "singleton": singleton,
        "total": total,
    }


def _log_fit(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    coeff = np.polyfit(np.log(np.clip(x, 1, None)), y, 1)
    return coeff[1] + coeff[0] * np.log(np.clip(x, 1, None))


def _stack(ax, strata: dict[str, np.ndarray], afr_start: int | None, ylabel: bool, arrows: bool = True) -> None:
    n = strata["total"].size
    x = np.arange(1, n + 1)
    stack = [strata[name] for name in FREQ_ORDER]
    ax.stackplot(
        x,
        *stack,
        colors=[FREQ_COLORS[name] for name in FREQ_ORDER],
        labels=[FREQ_LABELS[name] for name in FREQ_ORDER],
        lw=0,
        rasterized=True,
    )
    nonsing = strata["total"] - strata["singleton"]
    if afr_start is None:
        ax.plot(x, _log_fit(x, nonsing), color="black", lw=1.15, ls="-")
        ax.plot(x, _log_fit(x, strata["total"]), color="black", lw=1.15, ls="--")
    else:
        split = afr_start
        ax.plot(x[:split], _log_fit(x[:split], nonsing[:split]), color="black", lw=1.15, ls="-")
        ax.plot(x[split:], _log_fit(x[split:], nonsing[split:]), color="black", lw=1.15, ls="-")
        ax.plot(x[:split], _log_fit(x[:split], strata["total"][:split]), color="black", lw=1.15, ls="--")
        ax.plot(x[split:], _log_fit(x[split:], strata["total"][split:]), color="black", lw=1.15, ls="--")
        ax.axvline(split + 0.5, color="#666666", ls="--", lw=0.8, zorder=5)
        if arrows:
            ax.text(split * 0.5, 1.02, r"$\leftarrow$ Non-African", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.2)
            ax.text(split + (n - split) * 0.5, 1.02, r"African $\rightarrow$", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.2)
    ax.set_xlim(1, n)
    ax.set_ylim(0, strata["total"][-1] * 1.06)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: format_millions(v) if v >= 1000 else f"{v:.0f}"))
    if ylabel:
        ax.set_ylabel("Variant count (cumulative)")
    despine(ax)


def _ancestry_bar(ax, n: int) -> None:
    ax.set_xlim(1, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(left=False)
    for anc, start, end in _p2_blocks():
        ax.add_patch(
            Rectangle((start + 1, 0), end - start, 1.0, facecolor=ANC_COLORS[anc], edgecolor="none", zorder=2)
        )
        if end - start > 300:
            ax.text((start + end) / 2 + 0.5, 0.5, anc, ha="center", va="center", fontsize=5.6, color="white", fontweight="bold")
    ax.set_xlabel("Discovery sample (ancestry-ordered)")
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)


def _panel_yield(ax, new: np.ndarray) -> None:
    labels = []
    vals = []
    colors = []
    for anc, start, end in _p2_blocks():
        labels.append(anc)
        sl = slice(max(start, 1) if start == 0 else start, end)
        vals.append(new[sl].mean())
        colors.append(ANC_COLORS[anc])
    x = np.arange(len(labels))
    ax.bar(x, vals, color=colors, zorder=2, width=0.78)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("New unique SVs / sample")
    ax.set_title("Mean yield after ancestry-ordered addition", loc="left", pad=2)
    despine(ax)
    ax.text(0.02, 0.92, "First genome excluded from EAS mean", transform=ax.transAxes, fontsize=5.2, color="#666666")


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(21)
    blocks = _p2_blocks()
    afr_start = next(s for a, s, _ in blocks if a == "AFR")

    new_p2 = _incremental(P2_N, 26_000, float(SV_N["SV_50"]), afr_start, 1.70, rng)
    new_p1 = _incremental(P1_N, 22_000, 650_000.0, None, 1.0, rng)
    st_p2 = _strata(new_p2, shared_inf=1_200, major_plat=28_000, major_tau=80, sing_tau=1600, sing_max=0.50, afr_start=afr_start)
    st_p1 = _strata(new_p1, shared_inf=2_400, major_plat=18_000, major_tau=55, sing_tau=280, sing_max=0.48, afr_start=None)

    fig_p2 = new_figure(3.85)
    gs_p2 = GridSpec(2, 1, figure=fig_p2, left=0.10, right=0.98, top=0.88, bottom=0.12, hspace=0.08, height_ratios=[1.0, 0.10])
    ax_p2 = fig_p2.add_subplot(gs_p2[0])
    ax_bar_p2 = fig_p2.add_subplot(gs_p2[1], sharex=ax_p2)
    _stack(ax_p2, st_p2, afr_start, ylabel=True)
    ax_p2.set_title("Phase 2  ·  SVs ≥50 bp  ·  n = 12,261", loc="left", pad=8)
    ax_p2.tick_params(labelbottom=False)
    ax_p2.legend(loc="upper left", fontsize=5.5, bbox_to_anchor=(0.0, 0.98))
    n_ins = 80
    axins = ax_p2.inset_axes([0.34, 0.38, 0.30, 0.48])
    st_ins = {key: val[:n_ins] for key, val in st_p2.items()}
    _stack(axins, st_ins, None, ylabel=False, arrows=False)
    axins.set_title(f"First {n_ins} samples", loc="left", fontsize=5.8, pad=1)
    axins.tick_params(labelsize=5.0)
    axins.set_xlabel("")
    _ancestry_bar(ax_bar_p2, P2_N)
    save_panel(fig_p2, "figS8_phase2_discovery")

    def _p1(ax):
        _stack(ax, st_p1, None, ylabel=True)
        ax.set_title("Phase 1  ·  AFR  ·  n = 1,027", loc="left", pad=2)
        ax.set_xlabel("Discovery sample")

    ax_panel("figS8_phase1_discovery", 2.65, _p1, width=COL_IN * 1.15, left=0.16, right=0.96, top=0.86, bottom=0.16)
    ax_panel("figS8_yield", 2.65, lambda ax: _panel_yield(ax, new_p2), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)

    fig = new_figure(7.35)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.09,
        right=0.98,
        top=0.91,
        bottom=0.08,
        wspace=0.28,
        hspace=0.12,
        height_ratios=[2.15, 0.16, 1.35],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_bar = fig.add_subplot(gs[1, :], sharex=ax_a)
    ax_b = fig.add_subplot(gs[2, 0])
    ax_c = fig.add_subplot(gs[2, 1])

    _stack(ax_a, st_p2, afr_start, ylabel=True)
    ax_a.set_title("Phase 2  ·  SVs ≥50 bp  ·  n = 12,261", loc="left", pad=10)
    ax_a.tick_params(labelbottom=False)
    ax_a.legend(loc="upper left", fontsize=5.5, bbox_to_anchor=(0.0, 0.98))
    n_ins = 80
    axins = ax_a.inset_axes([0.34, 0.38, 0.30, 0.48])
    st_ins = {key: val[:n_ins] for key, val in st_p2.items()}
    _stack(axins, st_ins, None, ylabel=False, arrows=False)
    axins.set_title(f"First {n_ins} samples", loc="left", fontsize=5.8, pad=1)
    axins.tick_params(labelsize=5.0)
    axins.set_xlabel("")
    _ancestry_bar(ax_bar, P2_N)

    _stack(ax_b, st_p1, None, ylabel=True)
    ax_b.set_title("Phase 1  ·  AFR  ·  n = 1,027", loc="left", pad=2)
    ax_b.set_xlabel("Discovery sample")
    _panel_yield(ax_c, new_p2)

    stamp_mockup(fig, extra="Frequency classes recomputed after each added sample, as in Ebert et al. 2021; counts are placeholders")
    return save(fig, "figS8_ebert_discovery")
