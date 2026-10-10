"""Science-family figure style for AoU-LR Phase 2 mock-ups.

Print sizes follow Science 2-column figures (183 mm wide). Type is Helvetica
or Arial, at least MIN_FONT in the file so it stays at least 6 pt after the
draft manuscript scales the plate to 0.98 of the text width. Association
statistics are simulated unless a panel docstring says otherwise.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

# Science 2-column width (mm → inches). Draft TeX text width is 6.5 in
# (letter, 1 in margins); plates are included at 0.98\textwidth. Files are
# drawn at FIG_W and graphicx scales them down.
COL_IN = 3.54  # 90 mm single column
FIG_W = 7.20  # 183 mm double column
MAX_H = 9.40  # keep a figure on one letter page with caption
MANUSCRIPT_TEXT_IN = 6.5
PLACED_WIDTH_IN = 0.98 * MANUSCRIPT_TEXT_IN
PLACED_SCALE = PLACED_WIDTH_IN / FIG_W
# Science's floor for figure type in the printed size is 6 pt. MIN_FONT is
# the smallest size allowed in a figure file; 7 pt here is about 6.2 pt
# on the draft page.
MIN_PLACED_PT = 6.0
MIN_FONT = 7.0

# This package lives in the aou-lr-phase-2 code repo but reads its inputs from,
# and writes its figures into, the manuscript repo. Resolve that root
# explicitly: walking up from __file__ only worked while the code sat inside
# the manuscript checkout. Override with AOU_LR_MANUSCRIPT_ROOT when the two
# repos are not siblings.
REPO_ROOT = Path(__file__).resolve().parents[1]
_ENV_ROOT = os.environ.get("AOU_LR_MANUSCRIPT_ROOT")
ROOT = (
    Path(_ENV_ROOT).expanduser().resolve()
    if _ENV_ROOT
    else REPO_ROOT.parent / "aou-lr-phase-2-manuscript"
)
MANUSCRIPT_ROOT = ROOT
DATADIR = ROOT / "data"
OUTDIR = ROOT / "mockup_figures"
# Unlabeled one-plot files for mock-up remixing. Lettered combined plates stay
# in OUTDIR only for figures that already use real data.
PANELDIR = OUTDIR / "panels"
# Covariates ship with the code repo, so these are repo-relative.
COVARIATES_CANDIDATES = [
    REPO_ROOT / "covariates.v8.csv.gz",
    REPO_ROOT / "covariates.v7.csv.gz",
    REPO_ROOT / "scratch" / "covariates.v6.csv.gz",
]

# Okabe–Ito (colorblind-safe), plus a few neutrals.
OKABE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "pink": "#CC79A7",
    "grey": "#7A7A7A",
    "light": "#D9D9D9",
}

ANC_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "OTH"]
# Ancestry hues follow gnomAD GEN_ANC_COLORS, with saturation and lightness
# pulled in so the colours work as large fills rather than legend dots.
# MID is a quieter yellow-green than gnomAD's #33CC33 so it does not match EAS.
ANC_COLORS = {
    "AFR": "#8C4A9A",
    "AMR": "#C4473D",
    "EAS": "#1E8A5C",
    "EUR": "#4E8BB5",
    "MID": "#6E9B45",
    "SAS": "#E08A22",
    "OTH": "#8E9696",
}
# Same hues, darkened only where the fill fails as 6 pt text on white.
ANC_TEXT = {
    **ANC_COLORS,
    "EUR": "#44799D",
    "MID": "#597E38",
    "SAS": "#A46519",
    "OTH": "#6E7575",
}
ANC_LABELS = {
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "MID": "Middle Eastern",
    "SAS": "South Asian",
    "OTH": "Other / unassigned",
}

# Observed Phase 2 PacBio discovery counts (covariates.v6, n = 12,261).
P2_ANC_N = {
    "AFR": 2816,
    "OTH": 2419,
    "EUR": 2177,
    "AMR": 1966,
    "EAS": 1309,
    "SAS": 1197,
    "MID": 377,
}
P2_N = 12261
P1_N = 1027
P2_MID_N = 11128
P2_HIGH_N = 1133
P2_EHR_N = 10610
P2_RNA_N = 8327
P2_OLINK_N = 9243
P2_METH_N = 12261
P2_AIAN_N = 800
P2_ONT_N = 876

# Table 2 SV site counts (≥50 bp; ≥20 bp) from the current manuscript draft.
SV_N = {
    "DEL_50": 599_220,
    "DEL_20": 1_575_018,
    "INS_50": 1_875_567,
    "INS_20": 3_501_898,
    "INV_50": 41,
    "INV_20": 55,
    "SV_50": 2_474_828,
    "SV_20": 5_100_000,
}

SV_COLORS = {
    "DEL": OKABE["blue"],
    "INS": OKABE["orange"],
    "DUP": OKABE["green"],
    "INV": OKABE["pink"],
    "BND": OKABE["grey"],
}

MEI_COLORS = {
    "Alu": OKABE["sky"],
    "LINE-1": OKABE["vermillion"],
    "SVA": OKABE["green"],
}

# hg38 chromosome lengths (Mb) for Manhattan layouts.
HG38_MB = {
    1: 248.96,
    2: 242.19,
    3: 198.30,
    4: 190.21,
    5: 181.54,
    6: 170.81,
    7: 159.35,
    8: 145.14,
    9: 138.39,
    10: 133.80,
    11: 135.09,
    12: 133.28,
    13: 114.36,
    14: 107.04,
    15: 101.99,
    16: 90.34,
    17: 83.26,
    18: 80.37,
    19: 58.62,
    20: 64.44,
    21: 46.71,
    22: 50.82,
}


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": MIN_FONT,
            "axes.labelsize": MIN_FONT,
            "axes.titlesize": MIN_FONT + 0.5,
            "axes.titleweight": "regular",
            "axes.titlepad": 3.0,
            "xtick.labelsize": MIN_FONT,
            "ytick.labelsize": MIN_FONT,
            "legend.fontsize": MIN_FONT,
            "legend.frameon": False,
            "legend.handlelength": 1.2,
            "legend.handletextpad": 0.4,
            "legend.borderpad": 0.2,
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#222222",
            "axes.labelcolor": "#111111",
            "text.color": "#111111",
            "xtick.color": "#111111",
            "ytick.color": "#111111",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "lines.linewidth": 1.05,
            "lines.markersize": 3.5,
            "patch.linewidth": 0.5,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 600,
            "figure.dpi": 150,
            "mathtext.default": "regular",
        }
    )


def new_figure(height: float, width: float = FIG_W) -> plt.Figure:
    apply_style()
    fig = plt.figure(figsize=(width, min(height, MAX_H)))
    fig.patch.set_facecolor("white")
    return fig


def label_panel(ax: plt.Axes, letter: str, x: float = -0.02, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="right",
        va="bottom",
        color="#111111",
        clip_on=False,
    )


def stamp_mockup(fig: plt.Figure, extra: str = "") -> None:
    note = "MOCK-UP — placeholder statistics; not for submission"
    if extra:
        note = f"{note}. {extra}"
    fig.text(
        0.995,
        0.004,
        note,
        ha="right",
        va="bottom",
        fontsize=5.5,
        color="#8A8A8A",
        style="italic",
        transform=fig.transFigure,
    )


def despine(ax: plt.Axes, left: bool = True, bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


def thin_ticks(ax: plt.Axes, nbins: int = 4) -> None:
    ax.xaxis.set_major_locator(MaxNLocator(nbins))
    ax.yaxis.set_major_locator(MaxNLocator(nbins))


def save(fig: plt.Figure, stem: str, *, panel_stem: str | None = None) -> tuple[Path, Path]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pdf = OUTDIR / f"{stem}.pdf"
    png = OUTDIR / f"{stem}.png"
    fig.savefig(pdf, dpi=600, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.04)
    if panel_stem:
        save_panel(fig, panel_stem, stamp=False, close=False)
    plt.close(fig)
    return pdf, png


def save_panel(
    fig: plt.Figure,
    stem: str,
    *,
    stamp: bool = True,
    extra: str = "",
    close: bool = True,
) -> tuple[Path, Path]:
    """Write one unlabeled panel under mockup_figures/panels/."""
    if stamp:
        stamp_mockup(fig, extra=extra)
    PANELDIR.mkdir(parents=True, exist_ok=True)
    pdf = PANELDIR / f"{stem}.pdf"
    png = PANELDIR / f"{stem}.png"
    fig.savefig(pdf, dpi=600, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(png, dpi=400, bbox_inches="tight", pad_inches=0.05)
    if close:
        plt.close(fig)
    return pdf, png


def ax_panel(
    stem: str,
    height: float,
    drawer,
    *,
    width: float = FIG_W,
    left: float = 0.12,
    right: float = 0.97,
    top: float = 0.86,
    bottom: float = 0.14,
    projection: str | None = None,
    stamp: bool = True,
    extra: str = "",
) -> tuple[Path, Path]:
    """Draw one unlabeled panel on its own canvas."""
    fig = new_figure(height, width)
    ax = fig.add_axes([left, bottom, right - left, top - bottom], projection=projection)
    drawer(ax)
    return save_panel(fig, stem, stamp=stamp, extra=extra)


def subplot_panel(
    stem: str,
    height: float,
    drawer,
    *,
    width: float = FIG_W,
    left: float = 0.12,
    right: float = 0.97,
    top: float = 0.88,
    bottom: float = 0.12,
    stamp: bool = True,
    extra: str = "",
    **subplot_kw,
) -> tuple[Path, Path]:
    """Like ax_panel, but the host axes has a SubplotSpec (needed for nested grids)."""
    fig = new_figure(height, width)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    ax = fig.add_subplot(111, **subplot_kw)
    drawer(ax)
    return save_panel(fig, stem, stamp=stamp, extra=extra)


def covariates_path() -> Path | None:
    for path in COVARIATES_CANDIDATES:
        if path.is_file():
            return path
    return None


def chrom_offsets() -> dict[int, float]:
    """Cumulative Mb offsets with 8 Mb gaps for Manhattan x-axis."""
    off = {}
    x = 0.0
    for chrom, length in HG38_MB.items():
        off[chrom] = x
        x += length + 8.0
    return off


def chrom_centers() -> dict[int, float]:
    off = chrom_offsets()
    return {c: off[c] + HG38_MB[c] / 2.0 for c in HG38_MB}


def manhattan_x(chrom: np.ndarray, pos_mb: np.ndarray) -> np.ndarray:
    off = chrom_offsets()
    x = np.empty_like(pos_mb, dtype=float)
    for c in HG38_MB:
        mask = chrom == c
        x[mask] = off[c] + pos_mb[mask]
    return x


def simulate_gwas(
    rng: np.random.Generator,
    n_per_chrom: int = 900,
    peaks: list[tuple[int, float, float, float]] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return chrom, pos_Mb, -log10p with optional planted peaks.

    peaks: list of (chrom, pos_Mb, height, width_Mb)
    """
    chroms = []
    pos = []
    y = []
    for c, length in HG38_MB.items():
        n = max(80, int(n_per_chrom * length / 180.0))
        p = rng.uniform(0.05, length - 0.05, size=n)
        # Uniform p-values → exponential -log10p, clipped.
        u = rng.uniform(1e-3, 1.0, size=n)
        logp = -np.log10(u)
        logp = np.clip(logp, 0, 4.2)
        chroms.append(np.full(n, c))
        pos.append(p)
        y.append(logp)
    chrom = np.concatenate(chroms)
    pos_mb = np.concatenate(pos)
    logp = np.concatenate(y)
    if peaks:
        for c, mu, height, width in peaks:
            d = np.abs(pos_mb - mu)
            mask = (chrom == c) & (d < width * 3)
            bump = height * np.exp(-0.5 * ((pos_mb[mask] - mu) / width) ** 2)
            bump += rng.normal(0, 0.25, size=bump.size)
            logp[mask] = np.maximum(logp[mask], bump)
            # Guarantee a lead SNP near the mode.
            lead = np.argmin(np.abs(pos_mb - mu) + 1e3 * (chrom != c))
            logp[lead] = max(logp[lead], height + rng.uniform(-0.15, 0.25))
    return chrom, pos_mb, logp


def draw_manhattan(
    ax: plt.Axes,
    chrom: np.ndarray,
    pos_mb: np.ndarray,
    logp: np.ndarray,
    threshold: float | None = 7.3,
    s: float = 3.5,
) -> None:
    x = manhattan_x(chrom, pos_mb)
    for i, c in enumerate(HG38_MB):
        mask = chrom == c
        color = "#4C4C4C" if i % 2 == 0 else "#A6C8DC"
        ax.scatter(
            x[mask],
            logp[mask],
            s=s,
            c=color,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
    if threshold is not None:
        ax.axhline(threshold, color=OKABE["vermillion"], ls="--", lw=0.7, zorder=3)
    centers = chrom_centers()
    ticks = [centers[c] for c in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)]
    labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=5.5)
    ax.set_xlim(0, max(manhattan_x(np.array([22]), np.array([HG38_MB[22]]))) + 5)
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_xlabel("Chromosome")
    despine(ax)
    ax.margins(x=0.01)


def annotate_peaks(
    ax: plt.Axes,
    labels: list[tuple[int, float, float, str]],
) -> None:
    """Place locus names on a Manhattan plot.

    labels: (chrom, pos_Mb, y, text)
    """
    for chrom, pos_mb, y, text in labels:
        x = float(manhattan_x(np.array([chrom]), np.array([pos_mb]))[0])
        ax.annotate(
            text,
            xy=(x, y),
            xytext=(0, 4),
            textcoords="offset points",
            fontsize=5.6,
            ha="center",
            va="bottom",
            color="#111111",
            arrowprops=dict(arrowstyle="-", color="#444444", lw=0.45, shrinkA=0, shrinkB=1),
        )


def format_millions(n: float) -> str:
    if n >= 1_000_000:
        v = n / 1_000_000
        return f"{v:.1f}M" if v < 10 else f"{v:.0f}M"
    if n >= 1000:
        return f"{n / 1000:.0f}k"
    return f"{n:.0f}"
