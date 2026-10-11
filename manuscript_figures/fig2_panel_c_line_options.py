"""Figure 2 panel C: line alternatives to the stacked carrier-count bars.

Every option draws the same two series from data/fig2_panel_c_alts.json -- Phase 2
alleles a Phase 1 haplotype already carried, and alleles new in Phase 2 -- against
the number of haplotypes carrying the allele. Each option is drawn into an axes of
the size panel C actually gets in the plate (3.01 x 0.554 in), so the sheet shows
what the figure would show.

    python -m manuscript_figures.fig2_panel_c_line_options   (from the repo root)
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .fig2_options import MANUSCRIPT_ROOT
from .fig2_panel_c_options import ALTS_PATH, C_IN_P1, C_IN_P1_FILL, C_NEW, _CARRIER_LABELS, _carrier_counts
from .style import MIN_FONT, OUTDIR, apply_style, despine

FS = MIN_FONT
AX_W, AX_H = 3.01, 0.554  # panel C's real axes size in the plate
FINE_LABELS = ["1", "2", "3–4", "5–9", "10–19", "20–49", "50–99", "100–999", "1k–9,999", "≥10k"]
# Ticks that keep the 10-bin axis readable at 3 in wide.
FINE_TICKS = [0, 1, 3, 4, 5, 6, 7, 9]
LW = 1.15
MS = 2.6


def _fine() -> tuple[np.ndarray, np.ndarray]:
    d = json.loads(ALTS_PATH.read_text())
    return np.array(d["alleles_seen_in_phase1"], float), np.array(d["alleles_new_in_phase2"], float)


def _legend(ax, handles, **kw) -> None:
    opts = dict(
        frameon=False, fontsize=FS, handlelength=1.3, handletextpad=0.4, labelspacing=0.25,
        borderaxespad=0, loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, columnspacing=1.0,
    )
    ax.legend(handles=handles, **{**opts, **kw})


def _line_handles() -> list:
    kw = dict(lw=LW, marker="o", ms=MS, mec="none")
    return [
        Line2D([], [], color=C_IN_P1, label="In Phase 1", **kw),
        Line2D([], [], color=C_NEW, label="New in Phase 2", **kw),
    ]


def _fine_axis(ax, n: int) -> None:
    ax.set_xticks(FINE_TICKS)
    ax.set_xticklabels([FINE_LABELS[i] for i in FINE_TICKS])
    ax.set_xlim(-0.35, n - 0.65)
    ax.set_xlabel("Haplotypes carrying the allele")


def _plot_pair(ax, ys, yn) -> None:
    x = np.arange(ys.size)
    ax.plot(x, ys, color=C_IN_P1, lw=LW, marker="o", ms=MS, mec="none", zorder=3)
    ax.plot(x, yn, color=C_NEW, lw=LW, marker="o", ms=MS, mec="none", zorder=4)


# ---------------------------------------------------------------- options


def opt_bar_current(ax) -> None:
    """0 - the stacked bars now in the plate, for reference."""
    seen, new = _carrier_counts(json.loads(ALTS_PATH.read_text()))
    total = (seen + new).sum()
    ys, yn = 100.0 * seen / total, 100.0 * new / total
    x = np.arange(ys.size)
    ax.bar(x, ys, width=0.86, color=C_IN_P1_FILL, edgecolor=C_IN_P1, lw=0.6, zorder=2)
    ax.bar(x, yn, width=0.86, bottom=ys, color=C_NEW, edgecolor=C_NEW, lw=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(_CARRIER_LABELS)
    ax.set_xlim(-0.6, ys.size - 0.4)
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    from .fig2_panel_c_options import legend_handles

    _legend(ax, legend_handles(), handlelength=0.8)
    despine(ax)


def opt_lines_fine(ax) -> None:
    """1 - the same percentages unstacked, over the 10 bins the scan reports."""
    seen, new = _fine()
    total = (seen + new).sum()
    _plot_pair(ax, 100 * seen / total, 100 * new / total)
    _fine_axis(ax, seen.size)
    ax.set_ylim(0, None)
    ax.set_ylabel("Phase 2 alleles (%)")
    _legend(ax, _line_handles())
    despine(ax)


def opt_lines_coarse(ax) -> None:
    """2 - unstacked over the 6 bins the bars use now."""
    seen, new = _carrier_counts(json.loads(ALTS_PATH.read_text()))
    total = (seen + new).sum()
    _plot_pair(ax, 100 * seen / total, 100 * new / total)
    ax.set_xticks(np.arange(seen.size))
    ax.set_xticklabels(_CARRIER_LABELS)
    ax.set_xlim(-0.35, seen.size - 0.65)
    ax.set_ylim(0, None)
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax, _line_handles())
    despine(ax)


def opt_lines_within(ax) -> None:
    """3 - each series on its own denominator, so the two spectra are compared by shape."""
    seen, new = _fine()
    _plot_pair(ax, 100 * seen / seen.sum(), 100 * new / new.sum())
    _fine_axis(ax, seen.size)
    ax.set_ylim(0, None)
    ax.set_ylabel("Alleles in class (%)")
    _legend(ax, _line_handles())
    despine(ax)


def opt_lines_log(ax) -> None:
    """4 - allele counts on a log axis; the new series runs out, the known one does not."""
    seen, new = _fine()
    x = np.arange(seen.size)
    ax.plot(x, seen, color=C_IN_P1, lw=LW, marker="o", ms=MS, mec="none", zorder=3)
    ax.plot(x, np.maximum(new, 0.7), color=C_NEW, lw=LW, marker="o", ms=MS, mec="none", zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 3e7)
    ax.set_yticks([1, 1e3, 1e6])
    ax.set_yticklabels(["1", "10$^3$", "10$^6$"])
    _fine_axis(ax, seen.size)
    ax.set_ylabel("Alleles")
    _legend(ax, _line_handles())
    despine(ax)


def opt_share_new(ax) -> None:
    """5 - one line: the share of alleles at each carrier count that Phase 1 did not have."""
    seen, new = _fine()
    frac = 100 * new / (seen + new)
    x = np.arange(seen.size)
    ax.fill_between(x, 0, frac, color=C_NEW, alpha=0.18, lw=0, zorder=2)
    ax.plot(x, frac, color=C_NEW, lw=LW, marker="o", ms=MS, mec="none", zorder=4)
    ax.axhline(100 * new.sum() / (seen + new).sum(), color=C_IN_P1, lw=0.7, ls=(0, (3, 2)), zorder=3)
    ax.text(seen.size - 0.75, 100 * new.sum() / (seen + new).sum() + 4, "all alleles, 30%",
            fontsize=FS, color="#5A5A5A", ha="right", va="bottom")
    _fine_axis(ax, seen.size)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("New in Phase 2 (%)")
    despine(ax)


# Carrier counts spanned by each of the ten bins; the last runs to the call-rate floor.
_BIN_LO = np.array([1, 2, 3, 5, 10, 20, 50, 100, 1_000, 10_000], float)
_BIN_HI = np.array([1, 2, 4, 9, 19, 49, 99, 999, 9_999, 22_358], float)


def opt_density(ax) -> None:
    """6 - a true spectrum: alleles per carrier count, against carrier count, both log."""
    seen, new = _fine()
    width = _BIN_HI - _BIN_LO + 1
    mid = np.sqrt(_BIN_LO * _BIN_HI)
    ax.plot(mid, seen / width, color=C_IN_P1, lw=LW, marker="o", ms=MS, mec="none", zorder=3)
    ax.plot(mid, np.maximum(new / width, 1e-5), color=C_NEW, lw=LW, marker="o", ms=MS, mec="none", zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.7, 3e4)
    ax.set_ylim(1e-5, 3e7)
    ax.set_yticks([1e-4, 1, 1e4])
    ax.set_yticklabels(["10$^{-4}$", "1", "10$^4$"])
    ax.set_xticks([1, 10, 100, 1_000, 10_000])
    ax.set_xticklabels(["1", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("Haplotypes carrying the allele")
    ax.set_ylabel("Alleles per\ncarrier count")
    _legend(ax, _line_handles())
    despine(ax)


OPTIONS = (
    ("0 · stacked bars, as drawn now", opt_bar_current),
    ("1 · two lines, 10 bins, share of all Phase 2 alleles", opt_lines_fine),
    ("2 · two lines, the 6 bins the bars use", opt_lines_coarse),
    ("3 · two lines, each on its own denominator", opt_lines_within),
    ("4 · two lines, allele counts on a log axis", opt_lines_log),
    ("5 · one line, share new at each carrier count", opt_share_new),
    ("6 · two lines, alleles per carrier count, both axes log", opt_density),
)


def render(stem: str = "fig2c_line_options") -> Path:
    apply_style()
    left, right, top_pad, row_h = 0.72, 0.22, 0.46, 1.40
    fig_w = left + AX_W + right
    fig_h = len(OPTIONS) * row_h + 0.30
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    for i, (title, draw) in enumerate(OPTIONS):
        y_top = fig_h - 0.18 - i * row_h
        ax = fig.add_axes([left / fig_w, (y_top - top_pad - AX_H) / fig_h, AX_W / fig_w, AX_H / fig_h])
        draw(ax)
        fig.text(0.012, (y_top - 0.02) / fig_h, title, fontsize=7.5, fontweight="bold",
                 color="#111111", ha="left", va="top")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{stem}.png"
    fig.savefig(out, dpi=400)
    fig.savefig(OUTDIR / f"{stem}.pdf")
    plt.close(fig)
    return out


if __name__ == "__main__":
    print(render())
