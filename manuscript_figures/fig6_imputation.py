"""Figure 6: SV reference panel, imputation accuracy, and srWGS associations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    P1_N,
    P2_N,
    SV_N,
    annotate_peaks,
    ax_panel,
    despine,
    draw_manhattan,
    new_figure,
    save,
    simulate_gwas,
    stamp_mockup,
)


def _panel_panel_size(ax) -> None:
    labels = [
        "Discovery\nsamples",
        "Ancestry\ngroups",
        "SVs ≥50 bp\nin panel",
        "Frequency-\nfiltered SVs",
        "Imputed\nsrWGS N",
    ]
    p1 = np.array([P1_N, 1, 6.5e5, 1.8e5, 1.0e4], dtype=float)
    p2 = np.array([P2_N, 6, SV_N["SV_50"], 8.4e5, 2.45e5], dtype=float)
    # Plot as grouped bars on log scale for counts; ancestry groups linear — split.
    x = np.arange(len(labels))
    w = 0.36
    # Normalize each pair to Phase 2 = 1 for a fold-change view, annotate raw N.
    fold = p2 / p1
    ax.bar(x - w / 2, np.ones_like(x), width=w, color="#C8C8C8", label="Phase 1 (relative = 1)", zorder=2)
    ax.bar(x + w / 2, fold, width=w, color=OKABE["blue"], label="Phase 2 fold vs Phase 1", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fold change vs Phase 1")
    ax.set_title("Reference panel scale", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper left")
    raw_p2 = ["12,261", "6", "2.5M", "0.84M", "245k"]
    for i, (f, lab) in enumerate(zip(fold, raw_p2)):
        ax.text(i + w / 2, f + 0.4, lab, ha="center", va="bottom", fontsize=5.8, color="#222222")
    ax.set_ylim(0, max(fold) * 1.18)


def _panel_imputation_r2(ax, rng: np.random.Generator) -> None:
    maf = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50])
    # Phase 2 panel higher than Phase 1, especially for non-EUR rare SVs.
    def curve(base, slope, noise=0.015):
        y = 1 - np.exp(-slope * (maf * 100) ** 0.65) * (1 - base)
        y = np.clip(y + rng.normal(0, noise, maf.size), 0.05, 0.99)
        return y

    ax.plot(maf, curve(0.12, 0.55), color="#B0B0B0", lw=1.3, marker="o", ms=3.5, label="Phase 1 panel (AFR)")
    ax.plot(maf, curve(0.22, 0.85), color=OKABE["blue"], lw=1.4, marker="o", ms=3.5, label="Phase 2, AFR")
    ax.plot(maf, curve(0.18, 0.78), color=ANC_COLORS["AMR"], lw=1.3, marker="o", ms=3.5, label="Phase 2, AMR")
    ax.plot(maf, curve(0.16, 0.70), color=ANC_COLORS["EUR"], lw=1.3, marker="o", ms=3.5, label="Phase 2, EUR")
    ax.plot(maf, curve(0.10, 0.48), color=ANC_COLORS["EAS"], lw=1.3, marker="o", ms=3.5, label="Phase 2, EAS")
    ax.set_xscale("log")
    ax.set_xlim(0.0009, 0.55)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("SV minor allele frequency")
    ax.set_ylabel(r"Imputation dosage $r^{2}$")
    ax.set_title("Leave-one-out imputation accuracy (placeholder)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="lower right", fontsize=5.8)
    ax.axhline(0.8, color="#DDDDDD", ls=":", lw=0.6)


def _panel_enriched(ax) -> None:
    rows = [
        ("INS 312 bp  ·  chr8:11.6 Mb", "AFR", 0.082, 0.004, 9.4, "Neutrophil count"),
        ("DEL 1.8 kb  ·  HBA1/HBA2", "AFR", 0.041, 0.002, 11.2, "Mean corpuscular volume"),
        ("SVA VNTR  ·  chr6:28.9 Mb", "AMR", 0.063, 0.009, 8.1, "SLE phecode"),
        ("DEL 8.4 kb  ·  NPHS2 intron", "AFR", 0.012, 0.001, 7.9, "Urine ACR"),
        ("INS 6.0 kb L1  ·  chr19", "EAS", 0.037, 0.003, 7.4, "LDL cholesterol"),
        ("MEI Alu  ·  HLA-DRB1", "MID", 0.055, 0.011, 8.6, "Rheumatoid arthritis"),
        ("DUP 42 kb  ·  AMY1A", "SAS", 0.29, 0.18, 6.9, "Serum amylase"),
        ("DEL 3.1 kb  ·  GYPB", "AFR", 0.095, 0.008, 8.8, "MNS blood group / malaria"),
    ]
    y = np.arange(len(rows))
    ax.set_xlim(0, 10.2)
    ax.set_ylim(-0.55, len(rows) + 0.35)
    ax.axis("off")
    header_y = len(rows) + 0.05
    ax.text(0.05, header_y, "SV (placeholder)", fontsize=6, color="#666666")
    ax.text(5.15, header_y, "AF non-EUR / EUR", fontsize=6, color="#666666")
    ax.text(8.15, header_y, r"$-\log_{10}P$", fontsize=6, color="#666666")
    for i, (sv, anc, af_ne, af_eu, logp, trait) in enumerate(rows):
        yi = len(rows) - 1 - i
        ax.text(0.05, yi + 0.16, sv, fontsize=6.0, va="center", color="#222222")
        ax.text(0.05, yi - 0.16, f"{trait}  ·  {anc}", fontsize=5.5, va="center", color=ANC_COLORS[anc])
        ax.text(5.15, yi, f"{af_ne:.3f} / {af_eu:.3f}", fontsize=6.0, va="center")
        ax.barh(yi, logp / 12 * 1.7, left=8.15, height=0.42, color=ANC_COLORS[anc], zorder=2)
        ax.text(8.15 + logp / 12 * 1.7 + 0.08, yi, f"{logp:.1f}", fontsize=5.5, va="center")
    ax.set_title("Non-EUR–enriched imputed SV associations", loc="left", pad=6)


def _panel_phewas(ax, rng: np.random.Generator) -> None:
    chrom, pos, logp = simulate_gwas(
        rng,
        n_per_chrom=850,
        peaks=[
            (1, 196.7, 11.5, 0.45),  # CFH
            (6, 32.1, 16.8, 0.9),
            (6, 160.9, 13.1, 0.4),
            (8, 11.6, 9.2, 0.3),
            (16, 0.22, 12.0, 0.2),
            (17, 44.0, 10.4, 0.35),
            (19, 45.4, 9.8, 0.3),
            (22, 36.6, 8.9, 0.25),
        ],
    )
    draw_manhattan(ax, chrom, pos, logp, threshold=7.5, s=2.8)
    annotate_peaks(
        ax,
        [
            (6, 32.1, 17.0, "MHC"),
            (16, 0.22, 12.2, "HBA1/2"),
            (6, 160.9, 13.3, "LPA"),
            (1, 196.7, 11.7, "CFH"),
        ],
    )
    ax.set_ylim(0, 19.2)
    ax.set_title("Imputed SV-PheWAS in the AoU short-read cohort (placeholder)", loc="left", pad=2)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(11)
    ax_panel("fig6_panel_size", 2.85, _panel_panel_size, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("fig6_imputation_r2", 2.85, lambda ax: _panel_imputation_r2(ax, rng), width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.16)
    ax_panel("fig6_enriched", 2.55, _panel_enriched, left=0.08, right=0.98, top=0.88, bottom=0.14)
    ax_panel("fig6_imputed_phewas", 2.65, lambda ax: _panel_phewas(ax, rng), left=0.08, right=0.98, top=0.86, bottom=0.16)
    fig = new_figure(7.50)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.96,
        bottom=0.055,
        wspace=0.32,
        hspace=0.50,
        height_ratios=[1.0, 1.05, 1.15],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, :])

    _panel_panel_size(ax_a)
    _panel_imputation_r2(ax_b, rng)
    _panel_enriched(ax_c)
    _panel_phewas(ax_d, rng)

    stamp_mockup(fig)
    return save(fig, "fig6_imputation")
