"""Figure S5: Multi-omic QTL resource (eQTL, pQTL, meQTL, splicing)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.gridspec import GridSpec

from .style import (
    COL_IN,
    OKABE,
    annotate_peaks,
    ax_panel,
    despine,
    draw_manhattan,
    new_figure,
    save,
    simulate_gwas,
    stamp_mockup,
    subplot_panel,
)


def _panel_lead(ax) -> None:
    cats = ["eQTL\n(eGenes)", "pQTL\n(proteins)", "meQTL\n(CpG clusters)", "sQTL\n(introns)"]
    sv = np.array([0.31, 0.24, 0.18, 0.37])
    snv = 1.0 - sv
    x = np.arange(len(cats))
    ax.bar(x, sv, color=OKABE["orange"], label="SV is lead variant", zorder=2)
    ax.bar(x, snv, bottom=sv, color="#D5D5D5", label="SNV/indel is lead", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Fraction of QTLs")
    ax.set_title("SV as the lead causal variant (placeholder)", loc="left", pad=2)
    despine(ax)
    ax.legend(loc="upper right", fontsize=5.6)
    counts = ["[X] eGenes", "[X] proteins", "[X] CpGs", "[X] sQTLs"]
    for i, (v, lab) in enumerate(zip(sv, counts)):
        ax.text(i, v / 2, f"{v:.0%}", ha="center", va="center", fontsize=6.0, color="white")
        ax.text(i, 1.03, lab, ha="center", va="bottom", fontsize=5.4, color="#555555")


def _panel_coloc(ax) -> None:
    """Compact UpSet of eQTL ∩ pQTL ∩ meQTL ∩ EHR."""
    combos = {
        (1, 1, 1, 1): 42,
        (1, 1, 1, 0): 88,
        (1, 1, 0, 0): 95,
        (1, 0, 1, 0): 120,
        (1, 0, 0, 1): 61,
        (0, 1, 1, 0): 40,
        (1, 0, 0, 0): 310,
        (0, 0, 1, 0): 210,
    }
    keys = sorted(combos, key=combos.get, reverse=True)
    counts = np.array([combos[k] for k in keys], dtype=float)
    n = len(keys)
    fig = ax.figure
    spec = ax.get_subplotspec()
    ax.remove()
    inner = spec.subgridspec(2, 2, width_ratios=[0.28, 0.72], height_ratios=[0.62, 0.38], wspace=0.08, hspace=0.08)
    ax_bar = fig.add_subplot(inner[0, 1])
    ax_mat = fig.add_subplot(inner[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(inner[1, 0])

    ax_bar.bar(np.arange(n), counts, color="#3D3D3D", width=0.78, zorder=2)
    ax_bar.set_ylabel("Loci")
    ax_bar.tick_params(labelbottom=False)
    despine(ax_bar)
    ax_bar.set_title("Colocalizing SV signals across molecular and clinical layers", loc="left", pad=2)
    for i, c in enumerate(counts):
        ax_bar.text(i, c + 8, f"{int(c)}", ha="center", va="bottom", fontsize=5.2)

    labels = ["eQTL", "pQTL", "meQTL", "EHR"]
    mat = np.array(keys, dtype=float).T  # 4 x n
    ax_mat.set_xlim(-0.5, n - 0.5)
    ax_mat.set_ylim(-0.5, 3.5)
    ax_mat.set_yticks(range(4))
    ax_mat.set_yticklabels([])
    ax_mat.set_xticks([])
    for r in range(4):
        for c in range(n):
            on = mat[r, c] > 0
            ax_mat.scatter(c, r, s=22, c="#222222" if on else "#E2E2E2", zorder=3, linewidths=0)
        ax_mat.plot([-0.4, n - 0.6], [r, r], color="#EFEFEF", lw=5.5, zorder=1)
    for c, key in enumerate(keys):
        ys = [i for i, v in enumerate(key) if v]
        if len(ys) >= 2:
            ax_mat.plot([c, c], [min(ys), max(ys)], color="#222222", lw=0.85, zorder=2)
    despine(ax_mat, left=True, bottom=False)
    ax_mat.tick_params(bottom=False)

    set_n = [np.sum([combos[k] for k in keys if k[i]]) for i in range(4)]
    ax_set.barh(range(4), set_n, color="#4C4C4C", height=0.55)
    ax_set.set_yticks(range(4))
    ax_set.set_yticklabels(labels, fontsize=6.0)
    ax_set.yaxis.tick_right()
    ax_set.tick_params(axis="y", length=0, pad=3)
    ax_set.set_ylim(-0.5, 3.5)
    ax_set.invert_xaxis()
    despine(ax_set, left=False)
    ax_set.spines["right"].set_visible(False)
    ax_set.set_xticks([])
    return ax_bar


def _panel_splicing(ax, rng: np.random.Generator) -> None:
    n = 180
    length = rng.normal(24, 7.5, n)
    length = np.clip(length, 6, 52)
    psi = np.clip(0.18 + 0.014 * (length - 18) + rng.normal(0, 0.07, n), 0.02, 0.95)
    ax.scatter(length, psi, s=10, c=OKABE["green"], alpha=0.55, linewidths=0, rasterized=True)
    xs = np.linspace(8, 48, 40)
    ax.plot(xs, np.clip(0.18 + 0.014 * (xs - 18), 0.05, 0.9), color="#222222", lw=1.2)
    ax.set_xlabel("SVA VNTR copies at a placeholder locus")
    ax.set_ylabel("Exon inclusion (Ψ)")
    ax.set_title("Length-dependent splicing QTL (not insertion presence)", loc="left", pad=2)
    despine(ax)
    ax.text(0.05, 0.90, r"$r \approx 0.62$ (placeholder)", transform=ax.transAxes, fontsize=5.8, color="#444444")


def _panel_eqtl_manhattan(ax, rng: np.random.Generator) -> None:
    chrom, pos, logp = simulate_gwas(
        rng,
        n_per_chrom=520,
        peaks=[
            (6, 32.0, 18.2, 0.7),
            (19, 45.4, 12.4, 0.3),
            (1, 155.0, 10.1, 0.35),
            (16, 0.2, 9.4, 0.2),
            (11, 52.0, 8.8, 0.3),
        ],
    )
    draw_manhattan(ax, chrom, pos, logp, threshold=7.5, s=2.6)
    annotate_peaks(
        ax,
        [
            (6, 32.0, 18.4, "MHC eQTL"),
            (19, 45.4, 12.6, "APOE pQTL"),
            (16, 0.2, 9.6, "HBA meQTL"),
        ],
    )
    ax.set_ylim(0, 20.5)
    ax.set_title("SV-QTL Manhattan (placeholder lead SVs)", loc="left", pad=2)


def render() -> tuple[Path, Path]:
    rng = np.random.default_rng(17)
    ax_panel("figS5_lead", 2.85, _panel_lead, width=COL_IN * 1.12, left=0.16, right=0.96, top=0.88, bottom=0.18)
    subplot_panel("figS5_coloc", 3.15, _panel_coloc, width=COL_IN * 1.25, left=0.16, right=0.96, top=0.90, bottom=0.08)
    ax_panel("figS5_eqtl_manhattan", 2.55, lambda ax: _panel_eqtl_manhattan(ax, rng), left=0.08, right=0.98, top=0.86, bottom=0.16)
    ax_panel("figS5_splicing", 2.55, lambda ax: _panel_splicing(ax, rng), left=0.08, right=0.98, top=0.86, bottom=0.16)
    fig = new_figure(7.55)
    gs = GridSpec(
        3,
        2,
        figure=fig,
        left=0.08,
        right=0.98,
        top=0.95,
        bottom=0.07,
        wspace=0.32,
        hspace=0.50,
        height_ratios=[1.0, 1.15, 1.15],
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b_host = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, :])

    _panel_lead(ax_a)
    _panel_coloc(ax_b_host)
    _panel_eqtl_manhattan(ax_c, rng)
    _panel_splicing(ax_d, rng)

    stamp_mockup(fig)
    return save(fig, "figS5_qtl")
