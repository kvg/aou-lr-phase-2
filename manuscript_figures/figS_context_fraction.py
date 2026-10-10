"""Supplementary figure: context-class proportion versus length (Fig. 2A, read as shares).

The vertical axis is the share of sites in the bin, on a linear scale. Length bins
match Fig. 2A, with bins of fewer than 20 sites merged into a neighbor.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from .fig2_options import CTX, CTX_COLORS, FS, MIN_CELL, SPECS, load_region_bins
from .style import MANUSCRIPT_ROOT, despine, new_figure, save

BINS = MANUSCRIPT_ROOT / "data" / "fig2_length_spectrum.region.bins.tsv"


def _merge(edges: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Merge adjacent bins until each has at least MIN_CELL sites."""
    groups: list[tuple[int, int]] = []
    start = 0
    for j in range(len(counts)):
        if counts[start : j + 1].sum() >= MIN_CELL:
            groups.append((start, j + 1))
            start = j + 1
    if start < len(counts):
        if groups:
            groups[-1] = (groups[-1][0], len(counts))
        else:
            groups.append((0, len(counts)))
    new_edges = np.array([edges[a] for a, _b in groups] + [edges[groups[-1][1]]])
    merged = np.array([counts[a:b].sum(axis=0) for a, b in groups])
    return new_edges, merged


def _binned(region: dict, edges: np.ndarray, sign: int, partitions: tuple[str, ...]) -> np.ndarray:
    out = np.zeros((len(edges) - 1, len(CTX)))
    for (part, length), counts in region["counts"].items():
        if part not in partitions or length == 0 or np.sign(length) != sign:
            continue
        i = int(np.searchsorted(edges, abs(length), side="right")) - 1
        if 0 <= i < len(out):
            out[i] += counts
    return out


def render() -> tuple[Path, Path]:
    region = load_region_bins(BINS)
    fig = new_figure(3.15)
    gs = GridSpec(2, 3, figure=fig, left=0.07, right=0.985, top=0.84, bottom=0.16, wspace=0.08, hspace=0.22, width_ratios=[s["width"] for s in SPECS])
    handles = [Patch(facecolor=CTX_COLORS[k], label=k) for k in CTX]
    for row, (sign, name) in enumerate(((1, "INS"), (-1, "DEL"))):
        for col, spec in enumerate(SPECS):
            ax = fig.add_subplot(gs[row, col])
            parts = ("small", "sv") if spec["layer"] == "integrated" else ("ultralong",)
            counts = _binned(region, spec["edges"], sign, parts)
            edges, counts = _merge(spec["edges"], counts)
            total = counts.sum(axis=1)
            widths = np.diff(edges)
            bottom = np.zeros(len(total))
            frac = np.divide(counts, total[:, None], out=np.zeros_like(counts), where=total[:, None] >= MIN_CELL)
            keep = total >= MIN_CELL
            for j, key in enumerate(CTX):
                ax.bar(edges[:-1], np.where(keep, frac[:, j], 0.0), bottom=bottom, width=widths, align="edge", color=CTX_COLORS[key], lw=0, rasterized=True)
                bottom = bottom + np.where(keep, frac[:, j], 0.0)
            ax.set_xlim(*spec["xlim"])
            ax.set_ylim(0, 1)
            if spec["xscale"] == "log":
                ax.set_xscale("log")
            ax.set_xticks(spec["xticks"])
            ax.set_xticklabels(spec["labels"] if row == 1 else [])
            ax.set_yticks([0, 0.5, 1])
            ax.set_yticklabels(["0", "50", "100"] if col == 0 else [])
            for x, lab in spec.get("landmarks", ()):
                ax.axvline(x, color="#B8B8B8", lw=0.45, ls="--", zorder=1)
                if row == 0:
                    ax.text(x, 0.97, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=FS, color="#555555", clip_on=True)
            if row == 0:
                ax.text(0.5, 1.04, spec["title"], transform=ax.transAxes, ha="center", va="bottom", fontsize=FS)
            if col == 2:
                ax.text(0.98, 0.96, name, transform=ax.transAxes, ha="right", va="top", fontsize=FS, color="#444444", fontweight="bold")
            ax.tick_params(labelsize=FS)
            despine(ax)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, 1.0), ncol=4, frameon=False, fontsize=FS, handlelength=0.9, handleheight=0.7, columnspacing=0.8)
    fig.axes[0].set_ylabel("Sites (%)", fontsize=FS)
    fig.axes[3].set_ylabel("Sites (%)", fontsize=FS)
    fig.axes[4].set_xlabel("Variant length (bp)", fontsize=FS)
    pdf, png = save(fig, "figS_context_fraction")
    dest = MANUSCRIPT_ROOT / "figures"
    out = []
    for src in (pdf, png):
        target = dest / src.name
        shutil.copyfile(src, target)
        out.append(target)
    return out[0], out[1]


if __name__ == "__main__":
    print(render())
