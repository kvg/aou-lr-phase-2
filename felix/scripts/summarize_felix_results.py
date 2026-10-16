#!/usr/bin/env python3
"""Summarize FELIX step2 / FelixGenome merged result TSVs.

Used by FelixGenome.wdl (six-argument call, qq_cct.png / manhattan_cct.png outputs). FelixPilot.wdl's
Summarize task passes --p-column, --named-suffix and --report-title and expects the output set of
scripts/summarize_tractor_genome_results.py instead; see felix/README.md.

Handles FELIX admixed columns including:
  P_cct_admixed_c, P_hom_admixed_c, P_het_admixed_c, BETA_c_anc*, p.value_c_anc*, etc.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
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

        # chi2_1 = 2 * erfcinv(p)^2 is monotone DECREASING in p, so the median
        # of the chi-square values is the transform of the median p-value.
        # Only the one or two middle p-values need the scalar inverse, not all
        # n of them: the old loop took 49 s on 4M p-values and would have been
        # minutes per phenotype on the per-record arm, called twice.
        n = p.size
        if n == 0:
            return float("nan")
        mids = [p[n // 2]] if n % 2 else [p[n // 2 - 1], p[n // 2]]
        x = [2.0 * (erfcinv_one(float(v)) ** 2) for v in mids]
        return float(np.mean(x) / 0.4549364)

    # Genome-wide scatters of tens of millions of points are an opaque block
    # below the interesting range and are what drives this task's resident
    # memory (matplotlib's own allocation, not the dataframe). Keep every
    # point that could matter and thin the rest; lambda and the counts are
    # always computed on the full set, never on the thinned one.
    PLOT_MAX_POINTS = int(os.environ.get("FELIX_PLOT_MAX_POINTS", "300000"))
    PLOT_KEEP_BELOW_P = 1e-3

    def thin_for_plot(keep_all: np.ndarray, n: int) -> np.ndarray:
        """Boolean mask: everything in keep_all, plus a stride sample of the rest."""
        if n <= PLOT_MAX_POINTS:
            return np.ones(n, dtype=bool)
        mask = np.asarray(keep_all, dtype=bool).copy()
        budget = PLOT_MAX_POINTS - int(mask.sum())
        rest = np.flatnonzero(~mask)
        if budget > 0 and rest.size:
            # CEILING division. Floor made rest[::step] LONGER than the budget
            # (319,401 vs 298,411 on BI_164's repeat arm), and the [:budget]
            # that followed then dropped ~21,000 points from the END OF FILE
            # ORDER rather than sampling them away uniformly. ConcatPhenotypeScores
            # merges with LC_ALL=C sort -k1,1, which is byte-wise, so chr9 sorts
            # LAST -- the whole of chr9 and the top 46% of chr8 vanished from the
            # Manhattan dense band while the always-kept p<1e-3 points still drew
            # across the hole. Plot only: lambda and every count come from the
            # full vector, but the figure was wrong and looked like missing data.
            step = max(1, -(-rest.size // budget))
            mask[rest[::step][:budget]] = True
        return mask

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
            # exp/obs are computed on the FULL sorted vector so the quantiles
            # are right; only the plotted pairs are thinned.
            keep = thin_for_plot(p < PLOT_KEEP_BELOW_P, n)
            ax.scatter(exp[keep], obs[keep], s=4, alpha=0.45, c="#1f4e79",
                       rasterized=True)
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
            keep = thin_for_plot(
                d[pcol].to_numpy(dtype=float) < PLOT_KEEP_BELOW_P, len(d)
            )
            d = d.loc[keep]
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
        # Kept byte-for-byte equivalent to summarize_tractor_genome_results.
        # On Terra only one script is localized, so this fallback is ALWAYS the
        # one that runs; returning NaN here put NaN in calibration_summary.tsv
        # even though the WDL passes --pheno-cov.
        out = {
            "n_samples": np.nan,
            "n_cases": np.nan,
            "n_controls": np.nan,
            "n_missing_pheno": np.nan,
        }
        if pheno_cov is None or not Path(pheno_cov).exists():
            return out
        df = pd.read_csv(
            pheno_cov, sep="\t", usecols=lambda c: c in {"ID", phenotype}
        )
        if phenotype not in df.columns:
            return out
        y = pd.to_numeric(df[phenotype], errors="coerce")
        out["n_missing_pheno"] = int(y.isna().sum())
        y = y.dropna()
        out["n_samples"] = int(len(y))
        out["n_cases"] = int((y == 1).sum())
        out["n_controls"] = int((y == 0).sum())
        return out


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
    pheno_cov: Path | None,
    top_n: int,
    p_threshold: float,
) -> dict:
    safe = sanitize_pheno(phenotype)
    ph_dir = out_dir / "qc" / safe
    ph_dir.mkdir(parents=True, exist_ok=True)

    # No phenotype-named duplicate is written. It was a hardlink locally, which
    # cost nothing in the task, but it was ALSO a workflow output, so Cromwell
    # delocalized the same bytes a second time under a second name: 25.65 GB of
    # BI_164.felix.tsv beside an identical merged.felix.tsv, 24 TB of pure
    # duplication across 951 phenotypes on the per-record arm. Nothing in any
    # WDL, script or notebook ever consumed results_tsvs_named; the manifest
    # records the real path instead.

    # Column names from the header alone, then read only what each step needs.
    # This used to be a single pd.read_csv of the whole file, which is fine for
    # the repeat arms (~1 GB per phenotype) and fatal for the per-record arm
    # (~27 GB per phenotype against memory_gb = 16).
    head = pd.read_csv(path, sep="\t", nrows=0)
    pcol = pick_joint_p_column(head)
    chrom_col = "CHR" if "CHR" in head.columns else "CHROM"
    beta_cols = [c for c in head.columns if re.match(r"BETA_c_anc\d+$", c)]

    n_var = 0
    n_joint = 0
    lam_cct = float("nan")
    n_gw_sig = 0

    if pcol:
        # qq_plot needs the p-values; manhattan_genome needs chromosome,
        # position and the p-value. Nothing else is read over the full file.
        use = [c for c in (chrom_col, "POS", pcol) if c in head.columns]
        # No low_memory=False here. With usecols it is actively harmful: it
        # makes the C parser buffer the WHOLE file before dropping unwanted
        # columns, which is where this task's memory actually went -- 3.25 GB
        # to read three columns out of a 1.55 GB file, against 0.26 GB with
        # per-block parsing and declared dtypes. That is the difference
        # between surviving the per-record arm and not.
        dt = {}
        if chrom_col in use:
            dt[chrom_col] = "category"
        dt[pcol] = "float64"
        light = pd.read_csv(path, sep="\t", usecols=use, dtype=dt)
        n_var = len(light)
        pv = pd.to_numeric(light[pcol], errors="coerce")
        light[pcol] = pv
        n_joint = int(pv.notna().sum())
        lam_cct = qq_plot(
            pv.to_numpy(dtype=float), f"{safe} {pcol}", ph_dir / "qq_cct.png"
        )
        manhattan_genome(light, pcol, f"{safe} {pcol}", ph_dir / "manhattan_cct.png")
        n_gw_sig = int((pv.notna() & (pv < p_threshold)).sum())
        del light, pv
    else:
        # Result tables are gzipped by ConcatPhenotypeScores; pandas handles
        # that by extension but a bare open() would read the deflate stream.
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as fh:
            n_var = max(sum(1 for _ in fh) - 1, 0)

    # top_hits keeps every column, so stream the file and hold only top_n rows
    # rather than sorting the whole frame. Rows with a missing p-value are no
    # longer eligible; the old code could pad the table with them when fewer
    # than top_n markers were tested.
    if pcol:
        best = None
        for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=100_000):
            chunk[pcol] = pd.to_numeric(chunk[pcol], errors="coerce")
            chunk = chunk[chunk[pcol].notna()]
            if chunk.empty:
                continue
            best = chunk if best is None else pd.concat([best, chunk], ignore_index=True)
            best = best.nsmallest(top_n, pcol, keep="first")
        top_hits = best if best is not None else head.head(0)
    else:
        top_hits = pd.read_csv(path, sep="\t", nrows=top_n, low_memory=False)
    top_hits.to_csv(ph_dir / "top_hits.tsv", sep="\t", index=False)

    counts = phenotype_counts(pheno_cov, phenotype)
    return {
        "phenotype": phenotype,
        "safe_name": safe,
        "result_path": str(path),
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
    ap.add_argument("--arm", default=None,
                    help="arm label to record as the first column (main / ru / ru_biallelic). "
                         "Without it a calibration_summary.tsv lifted out of its tarball is "
                         "indistinguishable from the other arms' except by joint_p_col.")
    args = ap.parse_args()

    if len(args.results) != len(args.phenotypes):
        raise SystemExit("results and phenotypes length mismatch")

    out_dir = args.out_dir

    rows = []
    for path, pheno in zip(args.results, args.phenotypes):
        rows.append(
            summarize_one(
                path, pheno, out_dir, args.pheno_cov, args.top_n, args.p_threshold
            )
        )

    manifest = pd.DataFrame(rows)
    if args.arm:
        manifest.insert(0, "arm", args.arm)
    manifest.to_csv(out_dir / "results_manifest.tsv", sep="\t", index=False)

    cal_cols = [
        "phenotype",
        "joint_p_col",
        "lambda_gc_cct",
        "n_variants",
        "n_joint_tested",
        "n_gw_sig",
        "n_samples",
        "n_cases",
    ]
    if args.arm:
        cal_cols.insert(0, "arm")
    cal = manifest[cal_cols]
    cal.to_csv(out_dir / "calibration_summary.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'results_manifest.tsv'} ({len(rows)} phenotypes)")


if __name__ == "__main__":
    main()
