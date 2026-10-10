"""Figure 3: MHC/KIR haplotype graphs, frequencies, and disease associations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle
import matplotlib.pyplot as plt

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    ax_panel,
    despine,
    new_figure,
    save,
    stamp_mockup,
)

HLA_GENES = [
    ("HLA-A", 29.92, 0.38),
    ("HLA-C", 31.12, 0.22),
    ("HLA-B", 31.48, 0.28),
    ("DRB1", 32.52, 0.22),
    ("DQA1", 32.78, 0.18),
    ("DQB1", 32.99, 0.18),
    ("DPA1", 33.20, 0.16),
    ("DPB1", 33.38, 0.16),
]

# Allele colors are schematic and kept distinct from ancestry (Okabe–Ito).
ALLELE_PAL = [
    "#8DD3C7",
    "#BEBADA",
    "#FB8072",
    "#80B1D3",
    "#FDB462",
    "#B3DE69",
    "#FCCDE5",
    "#D9D9D9",
]

KIR_GENES = [
    "3DL3",
    "2DS2",
    "2DL2",
    "2DL3",
    "2DL5",
    "2DS3",
    "2DP1",
    "2DL1",
    "3DP1",
    "2DL4",
    "3DL1",
    "3DS1",
    "2DL5B",
    "2DS5",
    "2DS1",
    "2DS4",
    "3DL2",
]


def _mhc_haplotypes(rng: np.random.Generator):
    """Simulate 18 representative extended MHC haplotypes labeled by ancestry."""
    ancestries = (
        ["AFR"] * 5 + ["AMR"] * 3 + ["EUR"] * 4 + ["EAS"] * 3 + ["SAS"] * 2 + ["MID"] * 1
    )
    names = [
        "A*30-B*42-DRB1*03",
        "A*74-B*15-DRB1*13",
        "A*02-B*53-DRB1*11",
        "A*23-B*07-DRB1*15",
        "A*68-B*58-DRB1*12",
        "A*02-B*35-DRB1*08",
        "A*24-B*39-DRB1*04",
        "A*68-B*40-DRB1*14",
        "A*01-B*08-DRB1*03",
        "A*03-B*07-DRB1*15",
        "A*02-B*44-DRB1*04",
        "A*29-B*44-DRB1*07",
        "A*24-B*52-DRB1*15",
        "A*11-B*15-DRB1*12",
        "A*02-B*46-DRB1*09",
        "A*33-B*44-DRB1*07",
        "A*24-B*40-DRB1*15",
        "A*02-B*50-DRB1*07",
    ]
    # allele index per gene per haplotype
    alleles = rng.integers(0, 6, size=(len(names), len(HLA_GENES)))
    # Make some haplotypes share blocks (LD).
    alleles[8, :] = [0, 1, 0, 0, 0, 0, 2, 1]  # AH8.1
    alleles[9, :] = [2, 3, 2, 3, 1, 3, 0, 0]
    freqs = [
        0.018,
        0.014,
        0.022,
        0.011,
        0.009,
        0.031,
        0.019,
        0.016,
        0.068,
        0.041,
        0.037,
        0.022,
        0.044,
        0.028,
        0.019,
        0.025,
        0.017,
        0.021,
    ]
    return list(zip(names, ancestries, alleles, freqs))


def _panel_mhc(ax, rng: np.random.Generator) -> None:
    haps = _mhc_haplotypes(rng)
    n = len(haps)
    ax.set_xlim(29.82, 34.85)
    ax.set_ylim(-0.45, n + 1.7)
    for name, start, width in HLA_GENES:
        ax.add_patch(
            Rectangle((start, n + 0.35), width, 0.7, facecolor="#222222", edgecolor="none")
        )
        lift = 0.38 if name == "HLA-C" else 0.0
        ax.text(start + width / 2, n + 1.22 + lift, name, ha="center", va="bottom", fontsize=5.4)
    ax.text(34.15, n + 1.22, "AF", ha="center", va="bottom", fontsize=5.4, color="#555555")
    yticks = []
    ylabels = []
    for i, (name, anc, alleles, freq) in enumerate(haps):
        y = n - 1 - i
        yticks.append(y + 0.35)
        ylabels.append(name)
        ax.plot([29.85, 33.40], [y + 0.35, y + 0.35], color="#E6E6E6", lw=5.5, zorder=1)
        for j, (_, start, width) in enumerate(HLA_GENES):
            ax.add_patch(
                Rectangle(
                    (start, y),
                    width,
                    0.7,
                    facecolor=ALLELE_PAL[int(alleles[j]) % len(ALLELE_PAL)],
                    edgecolor="white",
                    linewidth=0.2,
                    zorder=2,
                )
            )
        ax.scatter(
            [33.58],
            [y + 0.35],
            s=16,
            c=ANC_COLORS[anc],
            zorder=4,
            edgecolors="white",
            linewidths=0.4,
        )
        bar_w = freq / 0.08 * 0.72
        ax.add_patch(
            Rectangle((33.72, y + 0.12), bar_w, 0.46, facecolor=ANC_COLORS[anc], edgecolor="none", zorder=2)
        )
        ax.text(33.72 + bar_w + 0.04, y + 0.35, f"{freq:.1%}", va="center", fontsize=5.0, color="#444444")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=5.2)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xlabel("chr6 position (Mb, schematic)")
    ax.set_title("Extended MHC haplotypes  ·  block color = allele (schematic)", loc="left", pad=2)
    despine(ax, left=False)
    ax.spines["left"].set_visible(False)
    for anc in ["AFR", "AMR", "EUR", "EAS", "SAS", "MID"]:
        ax.scatter([], [], s=16, c=ANC_COLORS[anc], edgecolors="white", linewidths=0.4, label=anc)
    ax.legend(loc="lower right", ncol=6, bbox_to_anchor=(1.0, 1.02), markerscale=1.0, fontsize=5.8, frameon=False)


def _panel_kir(ax, rng: np.random.Generator) -> None:
    n_hap = 14
    # KIR A (mostly inhibitory) vs B (activating) frameworks.
    presence = np.zeros((n_hap, len(KIR_GENES)), dtype=int)
    a_core = {"3DL3", "2DL3", "2DP1", "2DL1", "3DP1", "2DL4", "3DL1", "2DS4", "3DL2"}
    b_extra = {"2DS2", "2DL2", "2DL5", "2DS3", "3DS1", "2DL5B", "2DS5", "2DS1"}
    labels = []
    groups = []
    for i in range(n_hap):
        is_b = i >= 6
        groups.append("B" if is_b else "A")
        labels.append(f"A{i + 1}" if not is_b else f"B{i - 5}")
        for j, g in enumerate(KIR_GENES):
            if g in a_core:
                presence[i, j] = 1 if (not is_b or rng.random() > 0.15) else int(rng.random() > 0.4)
            elif g in b_extra:
                presence[i, j] = int(is_b and rng.random() > 0.25)
            else:
                presence[i, j] = int(rng.random() > 0.7)
    cmap = ListedColormap(["#F2F2F2", OKABE["green"]])
    ax.imshow(presence, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_xticks(np.arange(len(KIR_GENES)))
    ax.set_xticklabels(KIR_GENES, rotation=90, fontsize=5.5)
    ax.set_yticks(np.arange(n_hap))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title("KIR gene-content haplotypes", loc="left", pad=2)
    ax.axhline(5.5, color="#222222", lw=0.6)
    ax.text(-0.8, 2.5, "A", fontsize=8, fontweight="bold", va="center", color=OKABE["blue"])
    ax.text(-0.8, 10.0, "B", fontsize=8, fontweight="bold", va="center", color=OKABE["vermillion"])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
    ax.tick_params(length=0)
    ax.set_xlabel("KIR gene")


def _panel_miami(ax, rng: np.random.Generator) -> None:
    n = 180
    # HLA alleles (top) and KIR (bottom)
    hla_names = [
        "DQA1*05:01",
        "DQB1*02:01",
        "DRB1*03:01",
        "B*27:05",
        "C*06:02",
        "A*01:01",
        "MICA*007",
        "HLA-Y*02",
        "DPB1*04:01",
        "B*57:01",
        "C*07:01",
        "DRB1*15:01",
    ]
    kir_names = [
        "KIR3DL2*007",
        "KIR3DL1*001",
        "KIR2DS1*002",
        "KIR2DL1*003",
        "KIR3DS1*013",
        "KIR2DL3*001",
        "KIR2DS4*001",
        "KIR2DL2*001",
    ]
    x_hla = np.linspace(0.5, 11.5, n)
    y_hla = rng.exponential(0.6, n)
    y_hla = np.clip(y_hla, 0, 4)
    peaks_h = [1.2, 1.4, 3.1, 5.8, 8.9]
    heights_h = [18.5, 16.2, 12.4, 9.8, 8.1]
    for mu, h in zip(peaks_h, heights_h):
        y_hla = np.maximum(y_hla, h * np.exp(-0.5 * ((x_hla - mu) / 0.18) ** 2) + rng.normal(0, 0.2, n))
    x_kir = np.linspace(0.5, 11.5, n)
    y_kir = np.clip(rng.exponential(0.5, n), 0, 3.5)
    peaks_k = [2.0, 6.2, 9.4]
    heights_k = [11.2, 8.4, 7.6]
    for mu, h in zip(peaks_k, heights_k):
        y_kir = np.maximum(y_kir, h * np.exp(-0.5 * ((x_kir - mu) / 0.22) ** 2) + rng.normal(0, 0.2, n))

    ax.scatter(x_hla, y_hla, s=6, c=OKABE["blue"], alpha=0.7, linewidths=0, rasterized=True)
    ax.scatter(x_kir, -y_kir, s=6, c=OKABE["vermillion"], alpha=0.7, linewidths=0, rasterized=True)
    ax.axhline(0, color="#222222", lw=0.6)
    ax.axhline(7.3, color=OKABE["blue"], ls="--", lw=0.6)
    ax.axhline(-7.3, color=OKABE["vermillion"], ls="--", lw=0.6)
    ax.set_xlim(0, 12)
    ax.set_ylim(-13.8, 21.5)
    ax.set_ylabel("HLA (up)  /  KIR (down)")
    ax.set_xticks([1.2, 3.1, 5.8, 8.9])
    ax.set_xticklabels(["DQA1*05", "B*27:05", "C*06:02", "HLA-Y*02"], fontsize=5.5, rotation=25, ha="right")
    ax.set_xlabel("Ordered classical HLA alleles (KIR mirrored below)")
    ax.set_title("Miami plot of HLA and KIR associations", loc="left", pad=2)
    despine(ax)
    ax.text(1.2, 19.0, "T1D  DQA1/DQB1", fontsize=5.8, color=OKABE["blue"])
    ax.text(2.05, -12.4, "KIR3DL2*007  anemia", fontsize=5.8, color=OKABE["vermillion"])


def _panel_effects(ax, rng: np.random.Generator) -> None:
    traits = [
        "Type 1 diabetes",
        "Seroneg. spondylopathy",
        "Placental disorders",
        "Acute posthem. anemia",
        "Psoriasis",
        "Celiac disease",
        "Multiple sclerosis",
        "SLE",
    ]
    alleles = ["DQA1*05", "MICA*007", "HLA-Y*02", "KIR3DL2*007", "C*06:02", "DQB1*02", "DRB1*15", "DRB1*03"]
    ancs = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
    # z-like effect matrix
    z = rng.normal(0, 0.4, size=(len(traits), len(ancs)))
    planted = {
        (0, 3): 3.2,
        (0, 0): 2.1,
        (1, 3): 2.8,
        (2, 1): 2.4,
        (3, 0): 2.9,
        (4, 3): 2.6,
        (5, 3): 3.0,
        (6, 3): 2.2,
        (7, 0): 1.8,
    }
    for (i, j), v in planted.items():
        z[i, j] = v
    im = ax.imshow(z, cmap="coolwarm", vmin=-3.2, vmax=3.2, aspect="auto")
    ax.set_xticks(np.arange(len(ancs)))
    ax.set_xticklabels(ancs)
    ax.set_yticks(np.arange(len(traits)))
    ax.set_yticklabels([f"{t}  ·  {a}" for t, a in zip(traits, alleles)], fontsize=6)
    ax.set_title("Ancestry-stratified effect sizes (z)", loc="left", pad=2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
    ax.tick_params(length=0)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.ax.tick_params(labelsize=5.5)
    cbar.set_label("z", fontsize=6)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(3)
    ax_panel("fig3_mhc_haplotypes", 2.85, lambda ax: _panel_mhc(ax, rng), left=0.16, right=0.98, top=0.82, bottom=0.16)
    ax_panel("fig3_kir", 3.15, lambda ax: _panel_kir(ax, rng), width=COL_IN * 1.15, left=0.14, right=0.96, top=0.88, bottom=0.14)
    ax_panel("fig3_miami", 3.15, lambda ax: _panel_miami(ax, rng), width=COL_IN * 1.15, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("fig3_effects", 2.65, lambda ax: _panel_effects(ax, rng), left=0.16, right=0.92, top=0.88, bottom=0.16)
    fig = new_figure(7.70)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.18,
        right=0.98,
        top=0.93,
        bottom=0.06,
        wspace=0.28,
        hspace=0.50,
        height_ratios=[1.15, 1.05, 1.05],
        width_ratios=[1.25, 0.90],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])

    _panel_mhc(ax_a, rng)
    _panel_kir(ax_b, rng)
    _panel_miami(ax_c, rng)
    _panel_effects(ax_d, rng)

    stamp_mockup(fig)
    return save(fig, "fig3_mhc")
