#!/usr/bin/env python3
"""Summarize FELIX step2 / FelixGenome merged result TSVs.

Handles FELIX admixed columns including:
  P_cct_admixed_c, P_hom_admixed_c, P_het_admixed_c, BETA_c_anc*, p.value_c_anc*, etc.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Reuse helpers from tractor summarize when available.
try:
    from summarize_tractor_genome_results import (
        chrom_sort_key,
        lambda_gc_from_p,
        manhattan_genome,
        phenotype_counts,
        qq_plot,
        sanitize_pheno,
    )
except ImportError:
    def sanitize_pheno(name: str) -> str:
        s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
        return s or "phenotype"

    def chrom_sort_key(chrom: str) -> tuple:
        c = str(chrom).replace("chr", "").replace("CHR", "")
        if c.isdigit():
            return (0, int(c))
        order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        return (1, order.get(c.upper(), 100), c)

    def lambda_gc_from_p(pvals: np.ndarray) -> float:
        p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
        if p.size < 10:
            return float("nan")
        from math import erfc, sqrt

        def erfcinv_one(pi: float) -> float:
            lo, hi = 0.0, 10.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if erfc(mid) > pi:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)

        x = 2.0 * (np.array([erfcinv_one(float(v)) for v in p]) ** 2)
        return float(np.median(x) / 0.4549364)

    def qq_plot(pvals: np.ndarray, title: str, out_png: Path) -> float:
        p = np.sort(pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)])
        lam = lambda_gc_from_p(p)
        fig, ax = plt.subplots(figsize=(5, 5))
        if p.size == 0:
            ax.set_title(title + " (no p-values)")
        else:
            n = p.size
            exp = -np.log10(np.arange(1, n + 1) / (n + 1))
            obs = -np.log10(p)
            ax.scatter(exp, obs, s=4, alpha=0.45, c="#1f4e79", rasterized=True)
            m = max(float(np.nanmax(exp)), float(np.nanmax(obs)))
            ax.plot([0, m], [0, m], ls="--", c="gray", lw=1)
            ax.set_title(f"{title}\nλGC≈{lam:.3f}  n={n:,}")
        ax.set_xlabel(r"Expected $-\log_{10}(p)$")
        ax.set_ylabel(r"Observed $-\log_{10}(p)$")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return lam

    def manhattan_genome(df: pd.DataFrame, pcol: str, title: str, out_png: Path) -> None:
        d = df.dropna(subset=[pcol]).copy()
        d = d[np.isfinite(d[pcol]) & (d[pcol] > 0)]
        fig, ax = plt.subplots(figsize=(12, 3.5))
        if d.empty:
            ax.set_title(title + " (empty)")
        else:
            chrom_col = "CHR" if "CHR" in d.columns else "CHROM"
            d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
            d["nlp"] = -np.log10(d[pcol].astype(float))
            d = d.dropna(subset=["POS", "nlp"])
            chroms = sorted(d[chrom_col].astype(str).unique(), key=chrom_sort_key)
            offset = 0.0
            ticks = []
            for i, chrom in enumerate(chroms):
                sub = d.loc[d[chrom_col].astype(str) == chrom]
                x = sub["POS"].to_numpy(dtype=float) + offset
                ax.scatter(x, sub["nlp"], s=3, alpha=0.45, c="#1f4e79", rasterized=True)
                mid = offset + 0.5 * (sub["POS"].max() - sub["POS"].min() + 1)
                ticks.append((mid, chrom.replace("chr", "")))
                offset += float(sub["POS"].max()) + 1e6
            ax.set_xticks([t[0] for t in ticks])
            ax.set_xticklabels([t[1] for t in ticks], rotation=90, fontsize=7)
            ax.set_title(title)
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)

    def phenotype_counts(pheno_cov: Path | None, phenotype: str) -> dict:
        return {"n_samples": np.nan, "n_cases": np.nan, "n_controls": np.nan}


FELIX_JOINT_P_COLS = [
    "P_cct_admixed_c",
    "P_cct_admixed",
    "P_hom_admixed_c",
    "P_hom_admixed",
    "P_het_admixed_c",
    "P_het_admixed",
    "p.value_ancALL",
]


def pick_joint_p_column(df: pd.DataFrame) -> str | None:
    for c in FELIX_JOINT_P_COLS:
        if c in df.columns:
            return c
    for c in ["P", "pvalue", "p.value"]:
        if c in df.columns:
            return c
    return None


def summarize_one(
    path: Path,
    phenotype: str,
    out_dir: Path,
    named_dir: Path,
    pheno_cov: Path | None,
    top_n: int,
    p_threshold: float,
) -> dict:
    safe = sanitize_pheno(phenotype)
    ph_dir = out_dir / "qc" / safe
    ph_dir.mkdir(parents=True, exist_ok=True)

    named_path = named_dir / f"{safe}.felix.tsv"
    shutil.copy2(path, named_path)

    df = pd.read_csv(path, sep="\t", low_memory=False)
    pcol = pick_joint_p_column(df)
    n_var = len(df)
    n_joint = int(df[pcol].notna().sum()) if pcol else 0

    lam_cct = float("nan")
    n_gw_sig = 0
    if pcol:
        lam_cct = qq_plot(
            df[pcol].to_numpy(dtype=float),
            f"{safe} {pcol}",
            ph_dir / "qq_cct.png",
        )
        manhattan_genome(df, pcol, f"{safe} {pcol}", ph_dir / "manhattan_cct.png")
        n_gw_sig = int(
            (
                df[pcol].notna()
                & (pd.to_numeric(df[pcol], errors="coerce") < p_threshold)
            ).sum()
        )

    beta_cols = [c for c in df.columns if re.match(r"BETA_c_anc\d+$", c)]
    top = df.copy()
    if pcol:
        top = top.sort_values(pcol, na_position="last")
    top_hits = top.head(top_n)
    top_hits.to_csv(ph_dir / "top_hits.tsv", sep="\t", index=False)

    counts = phenotype_counts(pheno_cov, phenotype)
    return {
        "phenotype": phenotype,
        "safe_name": safe,
        "result_path": str(named_path),
        "n_variants": n_var,
        "n_joint_tested": n_joint,
        "joint_p_col": pcol or "",
        "lambda_gc_cct": lam_cct,
        "n_gw_sig": n_gw_sig,
        "n_beta_c_cols": len(beta_cols),
        **counts,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", nargs="+", required=True, type=Path)
    ap.add_argument("--phenotypes", nargs="+", required=True)
    ap.add_argument("--pheno-cov", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    args = ap.parse_args()

    if len(args.results) != len(args.phenotypes):
        raise SystemExit("results and phenotypes length mismatch")

    out_dir = args.out_dir
    named_dir = out_dir / "results_by_phenotype"
    named_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path, pheno in zip(args.results, args.phenotypes):
        rows.append(
            summarize_one(
                path, pheno, out_dir, named_dir, args.pheno_cov, args.top_n, args.p_threshold
            )
        )

    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "results_manifest.tsv", sep="\t", index=False)

    cal = manifest[
        [
            "phenotype",
            "joint_p_col",
            "lambda_gc_cct",
            "n_variants",
            "n_joint_tested",
            "n_gw_sig",
            "n_samples",
            "n_cases",
        ]
    ]
    cal.to_csv(out_dir / "calibration_summary.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'results_manifest.tsv'} ({len(rows)} phenotypes)")


if __name__ == "__main__":
    main()
