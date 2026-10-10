"""Figure 2: LD, SV–EHR associations, and multi-omic QTL localizations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.pyplot as plt

from .style import (
    ANC_COLORS,
    ANC_ORDER,
    COL_IN,
    OKABE,
    P2_N,
    annotate_peaks,
    ax_panel,
    despine,
    draw_manhattan,
    new_figure,
    save,
    simulate_gwas,
    stamp_mockup,
)


def _panel_overview(ax) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.2)
    ax.axis("off")
    boxes = [
        (0.15, 1.15, 1.9, 1.35, f"Phase 2 LRS\nn = {P2_N:,}", OKABE["blue"]),
        (2.45, 1.15, 2.15, 1.35, "Within-ancestry\nSAIGE mixed models", OKABE["green"]),
        (5.05, 1.15, 2.0, 1.35, "METAL\nfixed-effect meta", OKABE["orange"]),
        (7.5, 1.15, 2.25, 1.35, "Replication\nsrWGS · ONT · P1", OKABE["pink"]),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.04,rounding_size=0.08",
                facecolor=color,
                edgecolor="none",
                alpha=0.18,
                mutation_aspect=0.6,
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.04,rounding_size=0.08",
                facecolor="none",
                edgecolor=color,
                linewidth=0.9,
            )
        )
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=6.4, color="#222222")
    for x0, x1 in [(2.05, 2.45), (4.6, 5.05), (7.05, 7.5)]:
        ax.annotate(
            "",
            xy=(x1, 1.82),
            xytext=(x0, 1.82),
            arrowprops=dict(arrowstyle="-|>", color="#444444", lw=0.9, mutation_scale=8),
        )
    ax.text(
        0.15,
        0.35,
        "Covariates: age, sex, coverage, ancestry PCs  ·  Phenotypes: EHR phecodes, labs, medications",
        fontsize=6.2,
        color="#444444",
        ha="left",
        va="center",
    )
    ax.set_title("Association design", loc="left", pad=1)


def _panel_ld_scatter(ax, rng: np.random.Generator) -> None:
    n = 420
    # Narrative: at associated loci the lead SV often tags the causal haplotype
    # better than the lead SNV, so most points sit above the diagonal.
    r2_snv = rng.beta(4.2, 3.0, n)
    r2_sv = np.clip(r2_snv + rng.beta(2.2, 3.5, n) * 0.45 + rng.uniform(0.02, 0.12, n), 0, 0.99)
    ax.scatter(r2_snv, r2_sv, s=8, c="#4C4C4C", alpha=0.45, linewidths=0, rasterized=True)
    ax.plot([0, 1], [0, 1], ls="--", color="#888888", lw=0.7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Lead SNV $r^{2}$ to causal haplotype")
    ax.set_ylabel(r"Lead SV $r^{2}$ to same haplotype")
    ax.set_title("SV vs SNV tagging at associated loci", loc="left", pad=2)
    despine(ax)
    ax.text(0.05, 0.90, "Points above the diagonal:\nSV is the better tag", fontsize=5.8, color="#444444")


def _panel_forest(ax, rng: np.random.Generator) -> None:
    loci = [
        ("NPHS2 del", "AFR", -0.42, 0.08),
        ("HBA1/2 α-thal", "AFR", -0.31, 0.06),
        ("LPA KIV-2", "EUR", 0.55, 0.07),
        ("LPA KIV-2", "AFR", 0.22, 0.09),
        ("UGT1A1 promoter", "EAS", 0.38, 0.07),
        ("CFH 155 kb dup", "EUR", 0.29, 0.08),
        ("G6PD complex", "AFR", -0.47, 0.10),
        ("APOL1 G1/G2 tag SV", "AFR", 0.61, 0.11),
    ]
    y = np.arange(len(loci))
    for i, (name, anc, beta, se) in enumerate(loci):
        color = ANC_COLORS[anc]
        ax.plot([beta - 1.96 * se, beta + 1.96 * se], [i, i], color=color, lw=1.2)
        ax.scatter([beta], [i], s=18, color=color, zorder=3, edgecolors="white", linewidths=0.3)
    ax.axvline(0, color="#888888", lw=0.6, ls="--")
    ax.set_xlim(-0.85, 0.95)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{name}  ({anc})" for name, anc, *_ in loci], fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xlabel("Effect size (log-odds or s.d.)")
    ax.set_title("Ancestry-stratified effects (placeholder loci)", loc="left", pad=2)
    despine(ax)


def _panel_qtl_locus(ax, rng: np.random.Generator) -> None:
    x = np.linspace(0, 1.0, 220)  # Mb around a gene
    def locus_curve(mu, h, w, noise=0.12):
        y = h * np.exp(-0.5 * ((x - mu) / w) ** 2)
        y += rng.exponential(0.15, size=x.size) * 0.15
        return y

    e = locus_curve(0.42, 9.5, 0.05)
    p = locus_curve(0.42, 7.2, 0.06)
    m = locus_curve(0.39, 6.1, 0.08)
    ehr = locus_curve(0.42, 8.4, 0.045)
    ax.plot(x, e + 22, color=OKABE["blue"], lw=1.0, label="eQTL")
    ax.fill_between(x, 22, e + 22, color=OKABE["blue"], alpha=0.15, linewidth=0)
    ax.plot(x, p + 14, color=OKABE["green"], lw=1.0, label="pQTL")
    ax.fill_between(x, 14, p + 14, color=OKABE["green"], alpha=0.15, linewidth=0)
    ax.plot(x, m + 7, color=OKABE["orange"], lw=1.0, label="meQTL")
    ax.fill_between(x, 7, m + 7, color=OKABE["orange"], alpha=0.15, linewidth=0)
    ax.plot(x, ehr, color=OKABE["vermillion"], lw=1.0, label="EHR phecode")
    ax.fill_between(x, 0, ehr, color=OKABE["vermillion"], alpha=0.15, linewidth=0)
    ax.axvline(0.42, color="#222222", ls=":", lw=0.7)
    ax.scatter([0.42], [31.2], marker="v", s=28, color="#111111", zorder=4)
    ax.text(0.445, 31.2, "Lead SV", fontsize=5.8, color="#111111", va="center")
    ax.set_xlim(0, 1)
    ax.set_yticks([0, 7, 14, 22])
    ax.set_yticklabels(["EHR", "meQTL", "pQTL", "eQTL"])
    ax.set_xlabel("Position in 1 Mb window (placeholder locus)")
    ax.set_ylabel(r"$-\log_{10}(P)$  (offset tracks)")
    ax.set_title("Colocalizing SV across molecular and clinical traits", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", ncol=4, bbox_to_anchor=(1.0, 1.02))


def _panel_manhattan(ax, rng: np.random.Generator) -> None:
    chrom, pos, logp = simulate_gwas(
        rng,
        n_per_chrom=700,
        peaks=[
            (6, 32.0, 14.5, 0.8),
            (1, 155.0, 10.2, 0.4),
            (16, 0.2, 11.8, 0.25),
            (6, 161.0, 12.4, 0.5),
            (22, 36.6, 9.6, 0.3),
            (2, 234.0, 8.8, 0.35),
        ],
    )
    draw_manhattan(ax, chrom, pos, logp, threshold=7.5, s=3.2)
    annotate_peaks(
        ax,
        [
            (6, 32.0, 14.7, "MHC"),
            (16, 0.2, 12.0, "HBA1/2"),
            (6, 161.0, 12.6, "LPA"),
            (22, 36.6, 9.8, "APOL1"),
        ],
    )
    ax.set_ylim(0, 16.2)
    ax.set_title("SV–EHR associations (placeholder peaks at known complex loci)", loc="left", pad=2)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(7)
    ax_panel("fig2_design", 1.55, _panel_overview, left=0.03, right=0.99, top=0.88, bottom=0.06)
    ax_panel(
        "fig2_ld_scatter",
        3.15,
        lambda ax: _panel_ld_scatter(ax, rng),
        width=COL_IN * 1.08,
        left=0.18,
        right=0.96,
        top=0.88,
        bottom=0.16,
    )
    ax_panel(
        "fig2_forest",
        3.35,
        lambda ax: _panel_forest(ax, rng),
        width=COL_IN * 1.25,
        left=0.38,
        right=0.96,
        top=0.88,
        bottom=0.14,
    )
    ax_panel("fig2_manhattan", 2.55, lambda ax: _panel_manhattan(ax, rng), left=0.08, right=0.98, top=0.86, bottom=0.16)
    ax_panel("fig2_coloc", 2.55, lambda ax: _panel_qtl_locus(ax, rng), left=0.10, right=0.98, top=0.86, bottom=0.16)

    rng = np.random.default_rng(7)
    fig = new_figure(7.55)
    gs = GridSpec(
        4,
        2,
        figure=fig,
        left=0.12,
        right=0.985,
        top=0.97,
        bottom=0.05,
        wspace=0.34,
        hspace=0.55,
        height_ratios=[0.72, 1.15, 1.2, 1.25],
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])
    ax_e = fig.add_subplot(gs[3, :])

    _panel_overview(ax_a)
    _panel_ld_scatter(ax_b, rng)
    _panel_forest(ax_c, rng)
    _panel_manhattan(ax_d, rng)
    _panel_qtl_locus(ax_e, rng)

    stamp_mockup(fig)
    return save(fig, "fig2_associations")
