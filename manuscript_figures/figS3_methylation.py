"""Figure S3: Native HiFi methylation maps from pb-CpG-tools pileups."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .style import (
    OKABE,
    ROOT,
    despine,
    label_panel,
    new_figure,
    save,
)

MID_PASS_COV = 16.3
HIGH_PASS_COV = 33.2
CMAP = LinearSegmentedColormap.from_list("meth", ["#F4F8FB", OKABE["blue"]])

DATA_DIRS = (
    ROOT / "data" / "figS3",
    ROOT,
)


def _data_file(name: str) -> Path:
    for folder in DATA_DIRS:
        path = folder / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"Missing {name}. Expected under data/figS3/ or the manuscript root."
    )


def _load_pairs(path: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        next(handle)
        for raw in handle:
            parts = raw.split("\t")
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _load_stats(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        x_col = "coverage" if reader.fieldnames and "coverage" in reader.fieldnames else "combined_mean_cov"
        xs: list[float] = []
        ys: list[float] = []
        for row in reader:
            try:
                x = float(row[x_col])
                y = float(row["frac_hap_slots_ge5"])
            except (KeyError, TypeError, ValueError):
                continue
            xs.append(x)
            ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _load_concordance_json(path: Path) -> dict:
    return json.loads(path.read_text())


def _colorbar(fig, ax, mappable, label: str):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4.2%", pad=0.08)
    cb = fig.colorbar(mappable, cax=cax)
    cb.set_label(label, fontsize=6.0)
    cb.ax.tick_params(labelsize=5.5, width=0.5, length=2.0)
    cb.outline.set_linewidth(0.5)
    return cb


def _panel_concordance(fig, ax, xs: np.ndarray, ys: np.ndarray, meta: dict) -> None:
    hb = ax.hexbin(
        xs,
        ys,
        gridsize=70,
        bins="log",
        cmap=CMAP,
        mincnt=1,
        extent=(0, 100, 0, 100),
        linewidths=0,
        rasterized=True,
    )
    ax.plot([0, 100], [0, 100], ls="--", color="#888888", lw=0.7, zorder=3)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_box_aspect(1)
    ax.set_xlabel("Primrose site-mean 5mC (%)")
    ax.set_ylabel("Jasmine site-mean 5mC (%)")
    ax.set_title("Caller concordance (chr22 site means)", loc="left", pad=2)
    despine(ax)
    n = int(meta.get("n_sites") or xs.size)
    frac = float(meta.get("concordance_frac") or np.nan)
    r = float(meta.get("concordance_r") or np.nan)
    ax.text(
        0.04,
        0.96,
        rf"$r = {r:.3f}$  ·  {100 * frac:.1f}% with $|\Delta| < 10$ pp  ·  {n:,} sites",
        transform=ax.transAxes,
        fontsize=5.8,
        color="#333333",
        va="top",
    )
    _colorbar(fig, ax, hb, "CpGs")


def _panel_haplotype(fig, ax, cov: np.ndarray, frac: np.ndarray) -> None:
    hb = ax.hexbin(
        cov,
        frac,
        gridsize=50,
        bins="log",
        cmap=CMAP,
        mincnt=1,
        linewidths=0,
        rasterized=True,
    )
    median_y = float(np.median(frac))
    ax.axhline(median_y, color="#888888", ls=":", lw=0.7, zorder=3)
    ax.axvline(MID_PASS_COV, color=OKABE["blue"], ls=":", lw=0.7, zorder=3)
    ax.axvline(HIGH_PASS_COV, color=OKABE["orange"], ls=":", lw=0.7, zorder=3)
    ax.text(MID_PASS_COV + 0.6, 0.08, "mid-pass", fontsize=5.5, color=OKABE["blue"])
    ax.text(HIGH_PASS_COV + 0.6, 0.08, "high-pass", fontsize=5.5, color=OKABE["orange"])
    ax.set_xlabel("Mean HiFi coverage (×)")
    ax.set_ylabel("Haplotype-slot fraction ≥5×")
    ax.set_title("Haplotype-assigned 5mC versus coverage", loc="left", pad=2)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, max(58.0, float(np.nanpercentile(cov, 99.5)) + 2.0))
    despine(ax)
    ax.text(
        0.04,
        0.96,
        f"median {100 * median_y:.1f}%  ·  n = {cov.size:,}",
        transform=ax.transAxes,
        fontsize=5.8,
        color="#333333",
        va="top",
    )
    _colorbar(fig, ax, hb, "Samples")


def render() -> tuple[Path, Path]:
    xs, ys = _load_pairs(_data_file("methylation_concordance.pairs.tsv.gz"))
    cov, frac = _load_stats(_data_file("pbcpg_sample_stats.tsv"))
    meta = _load_concordance_json(_data_file("methylation_concordance.json"))
    fig = new_figure(3.55)
    gs = GridSpec(
        1,
        2,
        figure=fig,
        left=0.07,
        right=0.97,
        top=0.88,
        bottom=0.16,
        wspace=0.38,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    _panel_concordance(fig, ax_a, xs, ys, meta)
    label_panel(ax_a, "A", x=-0.14, y=1.08)
    _panel_haplotype(fig, ax_b, cov, frac)
    label_panel(ax_b, "B", x=-0.14, y=1.08)
    return save(fig, "figS3_methylation")
