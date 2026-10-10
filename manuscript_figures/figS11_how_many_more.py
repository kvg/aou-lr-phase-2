"""Figure S11: Is SV discovery plateauing, and how many sites remain?"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter

from .style import (
    COL_IN,
    OKABE,
    P2_N,
    SV_N,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    stamp_mockup,
)

# Placeholder ≥50 bp splits (same fractions as S9/S10).
# b: power-law yield exponent (larger → closer to plateau).
# f1, f2: singleton / doubleton fractions of the catalog (for Chao1).
# sat_frac: how far a saturating model would claim we are toward Smax.
CLASSES = (
    dict(name="INS unique", short="INS\nunique", n=0.22 * SV_N["INS_50"], b=0.72, f1=0.38, f2=0.16, sat=0.86, color=OKABE["blue"], ls="-"),
    dict(name="INS repetitive", short="INS\nrepetitive", n=0.78 * SV_N["INS_50"], b=0.34, f1=0.56, f2=0.12, sat=0.52, color=OKABE["orange"], ls="-"),
    dict(name="DEL unique", short="DEL\nunique", n=0.38 * SV_N["DEL_50"], b=0.78, f1=0.32, f2=0.18, sat=0.90, color=OKABE["sky"], ls="--"),
    dict(name="DEL repetitive", short="DEL\nrepetitive", n=0.62 * SV_N["DEL_50"], b=0.40, f1=0.48, f2=0.14, sat=0.58, color=OKABE["vermillion"], ls="--"),
)


def _yield_curve(cls: dict, n: np.ndarray) -> np.ndarray:
    """new unique sites at sample n under S ∝ n^(1-b)."""
    b = cls["b"]
    n_end = float(P2_N)
    # S(N) = c/(1-b) * (N**(1-b) - 1) = n_obs
    c = cls["n"] * (1.0 - b) / (n_end ** (1.0 - b) - 1.0)
    return c * np.power(n, -b)


def _chao1(n_obs: float, f1_frac: float, f2_frac: float) -> float:
    f1 = f1_frac * n_obs
    f2 = f2_frac * n_obs
    return n_obs + f1 * (f1 - 1.0) / (2.0 * (f2 + 1.0))


def _panel_yield(ax) -> None:
    n = np.unique(np.concatenate([np.linspace(20, 400, 40), np.geomspace(400, P2_N, 80)])).astype(float)
    for cls in CLASSES:
        y = _yield_curve(cls, n)
        ax.plot(n, y, color=cls["color"], ls=cls["ls"], lw=1.45, label=cls["name"])
        y_now = _yield_curve(cls, np.array([float(P2_N)]))[0]
        ax.scatter([P2_N], [y_now], s=16, color=cls["color"], zorder=4, linewidths=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(20, P2_N * 1.15)
    ax.set_xlabel("Discovery samples")
    ax.set_ylabel("New unique SVs / sample")
    ax.set_title("Yield is falling, but not to zero", loc="left", pad=2)
    despine(ax)
    ax.axvline(P2_N, color="#AAAAAA", ls=":", lw=0.7)
    ax.text(P2_N, 0.04, "n = 12,261", transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=5.5, color="#666666", rotation=90)
    ax.legend(loc="upper right", fontsize=5.5, ncol=2)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}" if v >= 1 else f"{v:.1f}"))


def _panel_remaining(ax) -> None:
    """Observed catalog, Chao1 unseen lower bound, saturating-model leftover."""
    x = np.arange(len(CLASSES))
    obs = np.array([c["n"] for c in CLASSES])
    chao_extra = np.array([_chao1(c["n"], c["f1"], c["f2"]) - c["n"] for c in CLASSES])
    sat_extra = np.array([c["n"] / c["sat"] - c["n"] for c in CLASSES])
    w = 0.36
    ax.bar(x - w / 2, obs / 1e6, width=w, color="#4A4A4A", label="Observed", zorder=2)
    ax.bar(x - w / 2, chao_extra / 1e6, width=w, bottom=obs / 1e6, color=OKABE["sky"], label="Chao1 unseen (lower bound)", zorder=2)
    ax.scatter(x + w / 2, (obs + sat_extra) / 1e6, s=28, marker="D", color=OKABE["vermillion"], zorder=4, label="Saturating-model total")
    ax.set_xticks(x)
    ax.set_xticklabels([c["short"] for c in CLASSES], fontsize=6.0)
    ax.set_ylabel("Sites (millions)")
    ax.set_title("How many more?  Bound, not a cap", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.3)
    ax.set_ylim(0, max((obs + np.maximum(chao_extra, sat_extra)) / 1e6) * 1.28)
    for i, (o, e) in enumerate(zip(obs, chao_extra)):
        ax.text(i - w / 2, (o + e) / 1e6 + 0.04, f"+{e / o:.0%}", ha="center", fontsize=5.5, color=OKABE["blue"])


def _panel_spectrum(ax) -> None:
    bins = ["1", "2", "3–5", "6–20", "21–100", ">100"]
    # Placeholder AC histograms (fractions). Repeats more singleton-heavy.
    specs = {
        "INS unique": np.array([0.38, 0.16, 0.18, 0.14, 0.09, 0.05]),
        "INS repetitive": np.array([0.56, 0.12, 0.13, 0.10, 0.06, 0.03]),
        "DEL unique": np.array([0.32, 0.18, 0.20, 0.15, 0.10, 0.05]),
        "DEL repetitive": np.array([0.48, 0.14, 0.16, 0.12, 0.07, 0.03]),
    }
    x = np.arange(len(bins))
    w = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5]
    for cls, off in zip(CLASSES, offsets):
        ax.bar(x + off * w, specs[cls["name"]], width=w, color=cls["color"], label=cls["name"], zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_xlabel("Allele count")
    ax.set_ylabel("Fraction of sites")
    ax.set_title("Why the bound is large: the catalog is still mostly rare", loc="left", pad=2)
    despine(ax)
    ax.set_ylim(0, 0.68)
    ax.legend(loc="upper right", fontsize=5.2, ncol=2)


def render() -> tuple[Path, Path]:
    ax_panel("figS11_yield", 2.55, _panel_yield, left=0.10, right=0.98, top=0.86, bottom=0.16)
    ax_panel("figS11_remaining", 2.85, _panel_remaining, width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("figS11_spectrum", 2.85, _panel_spectrum, width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)
    fig = new_figure(7.20)
    gs = GridSpec(
        2,
        2,
        figure=fig,
        left=0.09,
        right=0.98,
        top=0.93,
        bottom=0.09,
        wspace=0.28,
        hspace=0.42,
        height_ratios=[1.12, 1.0],
        width_ratios=[1.12, 1.0],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    _panel_yield(ax_a)
    _panel_remaining(ax_b)
    _panel_spectrum(ax_c)
    stamp_mockup(fig, extra="Chao1 and yield curves are layout placeholders from assumed singleton/doubleton fractions")
    return save(fig, "figS11_how_many_more")
