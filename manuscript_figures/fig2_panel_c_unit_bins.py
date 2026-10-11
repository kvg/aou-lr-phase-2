"""Figure 2 panel C with one bar per carrying-haplotype count.

Reads data/fig2_panel_c_carrier_hist.json, which holds the same two allele classes
as the binned scan -- Phase 2 alleles a Phase 1 haplotype already carried, and
alleles new in Phase 2 -- counted at one-haplotype resolution instead of in ten
bins. Every option is drawn into an axes the size panel C gets in the plate
(3.01 x 0.554 in).

    PYTHONPATH=. python3 -c "from manuscript_figures.fig2_panel_c_unit_bins import render; render()"
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from .fig2_options import MANUSCRIPT_ROOT
from .fig2_panel_c_options import C_IN_P1, C_IN_P1_FILL, C_NEW
from .style import MIN_FONT, OUTDIR, apply_style, despine

FS = MIN_FONT
AX_W, AX_H = 3.01, 0.554
HIST_PATH = MANUSCRIPT_ROOT / "data" / "fig2_panel_c_carrier_hist.json"
MAC_MIN = 20  # the count above which an allele can be tested on its own


def load() -> dict:
    d = json.loads(HIST_PATH.read_text())
    seen = np.array(d["seen_in_phase1"], float)
    new = np.array(d["new_in_phase2"], float)
    total = seen.sum() + new.sum() + d["seen_in_phase1_above_cap"] + d["new_in_phase2_above_cap"]
    return {"seen": seen, "new": new, "total": total, "raw": d}


def _legend(ax, **kw) -> None:
    handles = [
        Patch(facecolor=C_IN_P1_FILL, edgecolor=C_IN_P1, linewidth=0.5, label="In Phase 1"),
        Patch(facecolor=C_NEW, edgecolor=C_NEW, linewidth=0.5, label="New in Phase 2"),
    ]
    opts = dict(
        frameon=False, fontsize=FS, handlelength=0.8, handletextpad=0.4, labelspacing=0.25,
        borderaxespad=0, loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, columnspacing=1.0,
    )
    ax.legend(handles=handles, **{**opts, **kw})


def _stack(ax, x, ys, yn, width, lw=0.0) -> None:
    ax.bar(x, ys, width=width, color=C_IN_P1_FILL, edgecolor=C_IN_P1 if lw else "none", lw=lw, zorder=2)
    ax.bar(x, yn, width=width, bottom=ys, color=C_NEW, edgecolor=C_NEW if lw else "none", lw=lw, zorder=3)


def _cut(ax, n: int, y: float, label: str = f"{MAC_MIN} haplotypes") -> None:
    ax.axvline(MAC_MIN - 0.5, color="#B03A1A", lw=0.7, zorder=5)
    ax.text(MAC_MIN + 0.8, y, label, fontsize=FS, color="#B03A1A", ha="left", va="top")


# ---------------------------------------------------------------- options


def opt_unit_linear(ax, d: dict, last: int = 60) -> None:
    """1 - one bar per count to 60, the rest of the axis dropped."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    ys, yn = 100 * s / d["total"], 100 * n / d["total"]
    _stack(ax, np.arange(1, last + 1), ys, yn, 0.85)
    beyond = 100 * (d["total"] - d["seen"][1:last + 1].sum() - d["new"][1:last + 1].sum()) / d["total"]
    ax.set_xlim(0.3, last + 0.7)
    ax.set_ylim(0, None)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60])
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    ax.text(0.99, 0.95, f"{beyond:.0f}% of alleles are carried by > {last} haplotypes",
            transform=ax.transAxes, ha="right", va="top", fontsize=FS, color="#5A5A5A")
    _legend(ax)
    despine(ax)


def opt_unit_linear_tail(ax, d: dict, last: int = 60) -> None:
    """2 - the same, with everything above 60 gathered into a detached bar."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    ys, yn = 100 * s / d["total"], 100 * n / d["total"]
    x = np.arange(1, last + 1, dtype=float)
    _stack(ax, x, ys, yn, 0.85)
    ts = 100 * (d["seen"].sum() - s.sum() + d["raw"]["seen_in_phase1_above_cap"]) / d["total"]
    tn = 100 * (d["new"].sum() - n.sum() + d["raw"]["new_in_phase2_above_cap"]) / d["total"]
    xt = last + 4.0
    _stack(ax, np.array([xt]), np.array([ts]), np.array([tn]), 2.6)
    ax.set_xlim(0.3, xt + 2.2)
    ax.set_ylim(0, None)
    ax.set_xticks([1, 20, 40, 60, xt])
    ax.set_xticklabels(["1", "20", "40", "60", f"> {last}"])
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax)
    despine(ax)


def opt_unit_log_x(ax, d: dict) -> None:
    """3 - one bar per count across the whole range, on a log axis."""
    s, n = d["seen"][1:], d["new"][1:]
    x = np.arange(1, s.size + 1, dtype=float)
    ys, yn = 100 * s / d["total"], 100 * n / d["total"]
    _stack(ax, x, ys, yn, np.minimum(0.9, x * 0.9))
    ax.set_xscale("log")
    ax.set_xlim(0.55, 3e4)
    ax.set_ylim(0, None)
    ax.set_xticks([1, 10, 100, 1_000, 10_000])
    ax.set_xticklabels(["1", "10", "100", "1,000", "10,000"])
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax)
    despine(ax)


def opt_unit_log_x_zoom(ax, d: dict, last: int = 200) -> None:
    """4 - one bar per count to 200 on a log axis, where the bars stay separable."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    x = np.arange(1, last + 1, dtype=float)
    ys, yn = 100 * s / d["total"], 100 * n / d["total"]
    _stack(ax, x, ys, yn, np.minimum(0.9, x * 0.42))
    ax.set_xscale("log")
    ax.set_xlim(0.55, last * 1.25)
    ax.set_ylim(0, None)
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 200])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "50", "100", "200"])
    ax.set_ylabel("Phase 2 alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax)
    despine(ax)


def opt_composition(ax, d: dict, last: int = 60) -> None:
    """5 - one bar per count, each filled to 100%, so the bars show the split only."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    tot = np.maximum(s + n, 1)
    ys, yn = 100 * s / tot, 100 * n / tot
    _stack(ax, np.arange(1, last + 1), ys, yn, 0.85)
    ax.set_xlim(0.3, last + 0.7)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60])
    ax.set_ylabel("Alleles (%)")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax)
    despine(ax)


def opt_unit_log_y(ax, d: dict, last: int = 60) -> None:
    """6 - counts, not shares, on a log axis: the gray and slate fills overlap, not stack."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    x = np.arange(1, last + 1, dtype=float)
    ax.bar(x, s + n, width=0.85, color=C_IN_P1_FILL, edgecolor="none", zorder=2)
    ax.bar(x, np.maximum(n, 0.6), width=0.85, color=C_NEW, edgecolor="none", zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.6, 2e7)
    ax.set_yticks([1, 1e3, 1e6])
    ax.set_yticklabels(["1", "10$^3$", "10$^6$"])
    ax.set_xlim(0.3, last + 0.7)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60])
    ax.set_ylabel("Alleles")
    ax.set_xlabel("Haplotypes carrying the allele")
    _legend(ax)
    despine(ax)


OPTIONS = (
    ("1 · one bar per count to 60, axis stops there", opt_unit_linear),
    ("2 · one bar per count to 60, remainder as a detached bar", opt_unit_linear_tail),
    ("3 · one bar per count, whole range, log axis", opt_unit_log_x),
    ("4 · one bar per count to 200, log axis", opt_unit_log_x_zoom),
    ("5 · one bar per count, each filled to 100% (split only)", opt_composition),
    ("6 · counts on a log axis, fills overlap rather than stack", opt_unit_log_y),
)


def _mac_cut(ax, y: float = 0.5, va: str = "center") -> None:
    ax.axvline(MAC_MIN + 0.5, color="#B03A1A", lw=0.7, zorder=6)
    ax.text(MAC_MIN + 1.6, y, "testable alone\nat ≥ 20", fontsize=FS, color="#B03A1A",
            ha="left", va=va, linespacing=0.95, transform=ax.get_xaxis_transform())


def ref_composition(ax, d: dict, last: int = 60) -> None:
    """A - one bar per count, filled to 100%, with the 20-haplotype cut marked."""
    opt_composition(ax, d, last=last)
    _mac_cut(ax, y=0.30)


def ref_composition_100(ax, d: dict) -> None:
    """B - the same out to 100 haplotypes."""
    opt_composition(ax, d, last=100)
    ax.set_xticks([1, 20, 40, 60, 80, 100])
    _mac_cut(ax, y=0.30)


def ref_log_counts(ax, d: dict, last: int = 60) -> None:
    """C - allele counts on a log axis, gray total behind slate new."""
    opt_unit_log_y(ax, d, last=last)
    _mac_cut(ax, y=0.62)


def ref_composition_rate(ax, d: dict, last: int = 60) -> None:
    """D - composition bars, with how many alleles sit at each count as a line."""
    opt_composition(ax, d, last=last)
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    ax2 = ax.twinx()
    ax2.plot(np.arange(1, last + 1), s + n, color="#B03A1A", lw=0.9, zorder=6)
    ax2.set_yscale("log")
    ax2.set_ylim(1e4, 3e7)
    ax2.set_yticks([1e5, 1e6, 1e7])
    ax2.set_yticklabels(["10$^5$", "10$^6$", "10$^7$"])
    ax2.set_ylabel("Alleles", color="#B03A1A", rotation=270, labelpad=9)
    ax2.tick_params(axis="y", colors="#B03A1A")
    for side in ("top", "left"):
        ax2.spines[side].set_visible(False)
    ax2.spines["right"].set_color("#B03A1A")


REFINED = (
    ("A · unit bars filled to 100%, split only, to 60", ref_composition),
    ("B · the same, out to 100 haplotypes", ref_composition_100),
    ("C · unit bars, allele counts on a log axis", ref_log_counts),
    ("D · composition bars with an allele-count line", ref_composition_rate),
)


def _log_axis(ax, last: int, ticks: list[int]) -> None:
    ax.set_yscale("log")
    ax.set_ylim(300, 3e7)
    ax.set_yticks([1e3, 1e5, 1e7])
    ax.set_yticklabels(["10$^3$", "10$^5$", "10$^7$"])
    ax.set_xlim(0.3, last + 0.7)
    ax.set_xticks(ticks)
    ax.set_ylabel("Alleles")
    ax.set_xlabel("Haplotypes carrying the allele")
    despine(ax)


def _named_legend(ax, gray_label: str, **kw) -> None:
    handles = [
        Patch(facecolor=C_IN_P1_FILL, edgecolor=C_IN_P1, linewidth=0.5, label=gray_label),
        Patch(facecolor=C_NEW, edgecolor=C_NEW, linewidth=0.5, label="New in Phase 2"),
    ]
    opts = dict(
        frameon=False, fontsize=FS, handlelength=0.8, handletextpad=0.4, labelspacing=0.25,
        borderaxespad=0, loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, columnspacing=1.0,
    )
    ax.legend(handles=handles, **{**opts, **kw})


def fin_envelope(ax, d: dict, last: int = 60) -> None:
    """C1 - gray is every Phase 2 allele at that count, slate the new ones inside it."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    x = np.arange(1, last + 1, dtype=float)
    ax.bar(x, s + n, width=0.85, color=C_IN_P1_FILL, edgecolor="none", zorder=2)
    ax.bar(x, n, width=0.85, color=C_NEW, edgecolor="none", zorder=3)
    _log_axis(ax, last, [1, 10, 20, 30, 40, 50, 60])
    _named_legend(ax, "All Phase 2 alleles")


def fin_envelope_100(ax, d: dict) -> None:
    """C2 - the same out to 100 haplotypes."""
    fin_envelope(ax, d, last=100)
    ax.set_xticks([1, 20, 40, 60, 80, 100])


def fin_grouped(ax, d: dict, last: int = 60) -> None:
    """C3 - two bars per count, so gray keeps the 'In Phase 1' meaning it has in D."""
    s, n = d["seen"][1:last + 1], d["new"][1:last + 1]
    x = np.arange(1, last + 1, dtype=float)
    w = 0.42
    ax.bar(x - w / 2, s, width=w, color=C_IN_P1_FILL, edgecolor="none", zorder=2)
    ax.bar(x + w / 2, n, width=w, color=C_NEW, edgecolor="none", zorder=3)
    _log_axis(ax, last, [1, 10, 20, 30, 40, 50, 60])
    _named_legend(ax, "In Phase 1")


def fin_grouped_40(ax, d: dict) -> None:
    """C4 - the grouped bars at 40 counts, where each bar is twice as wide."""
    fin_grouped(ax, d, last=40)
    ax.set_xticks([1, 10, 20, 30, 40])


FINAL = (
    ("C1 · gray = all Phase 2 alleles, slate = the new ones, to 60", fin_envelope),
    ("C2 · the same, out to 100 haplotypes", fin_envelope_100),
    ("C3 · two bars per count, gray = in Phase 1, to 60", fin_grouped),
    ("C4 · two bars per count, to 40", fin_grouped_40),
)


def render(stem: str = "fig2c_unit_bin_options", options=None) -> Path:
    options = OPTIONS if options is None else options
    d = load()
    apply_style()
    left, right, top_pad, row_h = 0.72, 0.22, 0.46, 1.40
    fig_w = left + AX_W + right
    fig_h = len(options) * row_h + 0.30
    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    for i, (title, draw) in enumerate(options):
        y_top = fig_h - 0.18 - i * row_h
        ax = fig.add_axes([left / fig_w, (y_top - top_pad - AX_H) / fig_h, AX_W / fig_w, AX_H / fig_h])
        draw(ax, d)
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
