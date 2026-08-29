#!/usr/bin/env python3
"""Post-workflow helpers for TractorMixGenome QC / PheWAS-style plots.

Used by notebooks/terra/tractor_09_genome_post_workflow.ipynb. Not required inside the WDL.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_manifest(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_calibration(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def phewas_lambda_bar(summary: pd.DataFrame, out_png: Path, lam_col: str = "lambda_gc_joint_acpass") -> None:
    d = summary.sort_values(lam_col, ascending=False).copy()
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.35 * len(d))))
    ax.barh(d["phenotype"].astype(str), d[lam_col].astype(float), color="#4C78A8")
    ax.axvline(1.0, color="gray", ls="--", lw=1)
    ax.axvline(1.05, color="darkorange", ls=":", lw=1, label="1.05")
    ax.set_xlabel("λGC (AC-pass joint P)")
    ax.set_title("Tractor-Mix genomic inflation by phenotype")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def phewas_hit_counts(summary: pd.DataFrame, out_png: Path) -> None:
    d = summary.sort_values("n_genomewide_sig", ascending=False).copy()
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.35 * len(d))))
    ax.barh(
        d["phenotype"].astype(str),
        d["n_genomewide_sig"].astype(float),
        color="#F58518",
    )
    ax.set_xlabel("Genome-wide significant variants (AC-pass)")
    ax.set_title("Tractor-Mix hit counts by phenotype")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def phewas_locus_heatmap(
    hits: pd.DataFrame,
    out_png: Path,
    p_col: str = "P",
    max_loci: int = 40,
) -> None:
    """Phenotype × locus −log10(P) heatmap from phewas_genomewide_hits.tsv."""
    if hits.empty or p_col not in hits.columns:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.set_title("No genome-wide hits")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    d = hits.copy()
    chrom_col = "CHR" if "CHR" in d.columns else ("CHROM" if "CHROM" in d.columns else None)
    if chrom_col is None or "POS" not in d.columns:
        raise SystemExit("hits table needs CHR/CHROM and POS")
    d["locus"] = d[chrom_col].astype(str) + ":" + d["POS"].astype(str)
    if "ID" in d.columns:
        d["locus"] = d["locus"] + ":" + d["ID"].astype(str)
    d["nlp"] = -np.log10(pd.to_numeric(d[p_col], errors="coerce"))
    d = d.dropna(subset=["nlp", "phenotype"])

    # Keep top loci by min p across traits
    top_loci = (
        d.groupby("locus")["nlp"].max().sort_values(ascending=False).head(max_loci).index
    )
    d = d[d["locus"].isin(top_loci)]
    mat = d.pivot_table(index="phenotype", columns="locus", values="nlp", aggfunc="max")
    mat = mat.reindex(sorted(mat.index), axis=0)

    fig_w = max(8, 0.35 * mat.shape[1] + 2)
    fig_h = max(3.5, 0.4 * mat.shape[0] + 1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat.fillna(0).to_numpy(), aspect="auto", cmap="YlOrRd")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index.astype(str), fontsize=8)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns.astype(str), rotation=90, fontsize=6)
    ax.set_title(r"PheWAS-style heatmap (max $-\log_{10}P$ per locus)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label=r"$-\log_{10}(P)$")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compare_lambda_two_models(
    limited: pd.DataFrame,
    full: pd.DataFrame,
    out_png: Path,
    lam_col: str = "lambda_gc_joint_acpass",
) -> pd.DataFrame:
    a = limited.set_index("phenotype")[lam_col].rename("lambda_limited")
    b = full.set_index("phenotype")[lam_col].rename("lambda_full")
    merged = pd.concat([a, b], axis=1).dropna()
    merged["delta_full_minus_limited"] = merged["lambda_full"] - merged["lambda_limited"]

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(merged["lambda_limited"], merged["lambda_full"], s=40, c="#1f4e79")
    for pheno, row in merged.iterrows():
        ax.annotate(str(pheno), (row["lambda_limited"], row["lambda_full"]), fontsize=7)
    lims = [
        min(merged.min().min() * 0.98, 1.0),
        max(merged.max().max() * 1.02, 1.05),
    ]
    ax.plot(lims, lims, ls="--", c="gray", lw=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("λGC limited")
    ax.set_ylabel("λGC full")
    ax.set_title("Limited vs full covariate inflation")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return merged.reset_index()
