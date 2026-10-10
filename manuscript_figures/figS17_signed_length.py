"""Figure S17: Signed variant length — deletions left, SNVs center, insertions right.

Vertical lines at ±20 bp (DeepVariant / small-indel vs integrated SV) and
±10 kb (integrated SV vs ultralong). Lengths are placeholder draws scaled
to catalog sizes, as in Fig. S7.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.patches import Rectangle

from .figS7_size_spectrum import (
    C_DV,
    C_LARGE,
    C_MAIN,
    N_DV_DEL,
    N_DV_INS,
    N_LARGE,
    N_MAIN_DEL,
    N_MAIN_INS,
    N_SNV,
    _simulate,
)
from .style import (
    OKABE,
    ax_panel,
    despine,
    format_millions,
    new_figure,
    save,
    stamp_mockup,
)


def _signed_bins() -> np.ndarray:
    """1-bp bins through 20 bp, then log bins to 10 Mb, mirrored around 0."""
    inner = np.arange(1, 21, dtype=float)
    mid = np.logspace(np.log10(20.0), 4.0, 42)[1:]
    outer = np.logspace(4.0, 7.0, 30)[1:]
    pos = np.concatenate([inner, mid, outer])
    return np.concatenate([-pos[::-1], np.array([-0.45, 0.45]), pos])


def _hist(values: np.ndarray, n_true: int, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    raw, edges = np.histogram(values, bins=bins)
    scale = n_true / max(values.size, 1)
    return raw * scale, edges


def _region_rail(ax) -> None:
    ax.set_xlim(-3.2e6, 3.2e6)
    ax.set_xscale("symlog", linthresh=20, linscale=0.55)
    ax.set_ylim(0, 1)
    ax.axis("off")
    bands = [
        (-20, 20, C_DV, "SNV / small indel  ·  DeepVariant  ·  0–20 bp"),
        (-10_000, -20, C_MAIN, "Integrated SV"),
        (20, 10_000, C_MAIN, "Integrated SV  ·  20 bp – 10 kb"),
        (-3.2e6, -10_000, C_LARGE, "Ultralong"),
        (10_000, 3.2e6, C_LARGE, "Ultralong  ·  >10 kb"),
    ]
    for lo, hi, color, _lab in bands:
        ax.add_patch(Rectangle((lo, 0.18), hi - lo, 0.62, facecolor=color, edgecolor="none", alpha=0.92, zorder=2))
    ax.text(0, 0.49, "0–20 bp", ha="center", va="center", fontsize=6.0, color="white", fontweight="bold", zorder=3)
    ax.text(-400, 0.49, "SV", ha="center", va="center", fontsize=6.0, color="white", fontweight="bold", zorder=3)
    ax.text(400, 0.49, "SV", ha="center", va="center", fontsize=6.0, color="white", fontweight="bold", zorder=3)
    ax.text(-1.2e5, 0.49, ">10 kb", ha="center", va="center", fontsize=6.0, color="white", fontweight="bold", zorder=3)
    ax.text(1.2e5, 0.49, ">10 kb", ha="center", va="center", fontsize=6.0, color="white", fontweight="bold", zorder=3)
    ax.set_title("Callset partitions on a signed length axis", loc="left", pad=2)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(9)
    data = _simulate(rng)
    bins = _signed_bins()

    h_del, edges = _hist(-data["dv_del"], N_DV_DEL, bins)
    h_del2, _ = _hist(-data["main_del"], N_MAIN_DEL, bins)
    h_del3, _ = _hist(-data["large_del"], int(N_LARGE * 0.62), bins)
    h_ins, _ = _hist(data["dv_ins"], N_DV_INS, bins)
    h_ins2, _ = _hist(data["main_ins"], N_MAIN_INS, bins)
    h_ins3, _ = _hist(data["large_ins"], int(N_LARGE * 0.18), bins)
    h_left = h_del + h_del2 + h_del3
    h_right = h_ins + h_ins2 + h_ins3

    # SNVs occupy the bin that straddles 0.
    centers = 0.5 * (edges[:-1] + edges[1:])
    zero = int(np.argmin(np.abs(centers)))
    h_snv = np.zeros_like(h_left)
    h_snv[zero] = N_SNV

    def _hist_panel(ax):
        ax.axvspan(-20, 20, color=C_DV, alpha=0.10, zorder=0)
        ax.axvspan(-10_000, -20, color=C_MAIN, alpha=0.08, zorder=0)
        ax.axvspan(20, 10_000, color=C_MAIN, alpha=0.08, zorder=0)
        ax.axvspan(-3.2e6, -10_000, color=C_LARGE, alpha=0.10, zorder=0)
        ax.axvspan(10_000, 3.2e6, color=C_LARGE, alpha=0.10, zorder=0)
        ax.stairs(h_left / 1e6, edges, fill=True, color=OKABE["blue"], alpha=0.88, zorder=2, label="Deletions")
        ax.stairs(h_right / 1e6, edges, fill=True, color=OKABE["orange"], alpha=0.88, zorder=2, label="Insertions")
        ax.stairs(h_snv / 1e6, edges, fill=True, color="#4A4A4A", alpha=0.95, zorder=3, label="SNVs")
        for x in (-10_000, -20, 20, 10_000):
            ax.axvline(x, color="#222222", lw=0.7, ls="--", zorder=4)
        ax.axvline(0, color="#222222", lw=0.55, zorder=4)
        ax.set_xscale("symlog", linthresh=20, linscale=0.55)
        ax.set_xlim(-3.2e6, 3.2e6)
        ax.set_yscale("log")
        ax.set_ylim(5e-4, 2.2e2)
        ticks = [-1e6, -1e4, -20, 0, 20, 1e4, 1e6]
        ax.set_xticks(ticks)
        ax.set_xticklabels(["-1 Mb", "-10 kb", "-20 bp", "0", "20 bp", "10 kb", "1 Mb"])
        ax.set_ylabel("Sites (millions / bin)")
        ax.set_xlabel("Signed variant length  (deletions left, insertions right)")
        ax.set_title("Length spectrum across SNV, SV, and ultralong callsets", loc="left", pad=2)
        despine(ax)
        ax.legend(loc="upper left", fontsize=6.0)
        ax.text(-8, 80, format_millions(N_SNV), ha="center", va="bottom", fontsize=5.8, color="#222222")
        ax.annotate(
            "Alu",
            xy=(310, 0.08),
            xytext=(310, 0.8),
            fontsize=5.4,
            ha="center",
            color="#555555",
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )
        ax.annotate(
            "Alu",
            xy=(-310, 0.08),
            xytext=(-310, 0.8),
            fontsize=5.4,
            ha="center",
            color="#555555",
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )

    ax_panel("figS17_partition_rail", 1.15, _region_rail, left=0.10, right=0.98, top=0.78, bottom=0.18)
    ax_panel("figS17_signed_spectrum", 3.35, _hist_panel, left=0.10, right=0.98, top=0.88, bottom=0.16)

    fig = new_figure(5.15)
    ax_r = fig.add_axes([0.10, 0.82, 0.88, 0.12])
    ax = fig.add_axes([0.10, 0.12, 0.88, 0.66])
    _region_rail(ax_r)
    _hist_panel(ax)

    stamp_mockup(fig)
    return save(fig, "figS17_signed_length")
