#!/usr/bin/env python3
"""Summarize Tractor-Mix or FELIX association result TSVs.

Writes phenotype-named association tables, per-trait QQ / Manhattan / top hits,
a results manifest, and cohort-level λGC tables suitable for Terra outputs and
post-workflow notebooks.

Expected Tractor-Mix columns include:
  CHR POS ID REF ALT Chi2 P Eff_anc* SE_anc* Pval_anc* AC_count* include_anc*

FELIX Step 2 columns (pass --p-column P_cct_admixed_c --named-suffix .felix.tsv):
  CHR POS MarkerID Allele1 Allele2 P_cct_admixed_c P_het_admixed_c
  P_hom_admixed_c p.value_c_anc* BETA_c_anc* SE_c_anc*
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def erfcinv(p: np.ndarray) -> np.ndarray:
    try:
        from scipy.special import erfcinv as s_erfcinv

        return np.asarray(s_erfcinv(p), dtype=float)
    except Exception:
        out = np.empty(len(p), dtype=float)
        for i, pi in enumerate(np.asarray(p, dtype=float)):
            lo, hi = 0.0, 10.0
            target = float(pi)
            if target <= 0:
                out[i] = hi
                continue
            if target >= 2:
                out[i] = 0.0
                continue
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if math.erfc(mid) > target:
                    lo = mid
                else:
                    hi = mid
            out[i] = 0.5 * (lo + hi)
        return out


def lambda_gc_from_p(pvals: np.ndarray) -> float:
    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size < 10:
        return float("nan")
    x = 2.0 * (erfcinv(p) ** 2)
    return float(np.median(x) / 0.4549364)


def sanitize_pheno(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return s or "phenotype"


def chrom_sort_key(chrom: str) -> tuple:
    c = str(chrom).replace("chr", "").replace("CHR", "")
    if c.isdigit():
        return (0, int(c))
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (1, order.get(c.upper(), 100), c)


def ac_pass_mask(df: pd.DataFrame, pcol: str = "P") -> pd.Series:
    include_cols = [c for c in df.columns if c.startswith("include_anc")]
    if include_cols:
        inc = df[include_cols].astype(str).apply(
            lambda s: s.str.lower().isin(["true", "1", "1.0", "t"])
        )
        return inc.any(axis=1)
    if pcol in df.columns:
        return df[pcol].notna()
    if "P" in df.columns:
        return df["P"].notna()
    return pd.Series(False, index=df.index)


def default_extra_p_columns(df: pd.DataFrame, pcol: str) -> list[str]:
    """QQ extras for FELIX _c_ columns without changing Tractor-Mix defaults."""
    extras: list[str] = []
    if pcol == "P_cct_admixed_c":
        for c in ("P_het_admixed_c", "P_hom_admixed_c"):
            if c in df.columns:
                extras.append(c)
        extras.extend(
            sorted(
                c
                for c in df.columns
                if c.startswith("p.value_c_anc") and c not in extras
            )
        )
    return extras


def qq_plot(pvals: np.ndarray, title: str, out_png: Path) -> float:
    p = np.sort(pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)])
    lam = lambda_gc_from_p(p)
    fig, ax = plt.subplots(figsize=(5, 5))
    if p.size == 0:
        ax.set_title(title + " (no p-values)")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return lam
    n = p.size
    # Thin dense null for genome-wide plots
    if n > 200_000:
        idx = np.unique(
            np.concatenate(
                [
                    np.linspace(0, n - 1, 150_000, dtype=int),
                    np.arange(max(0, n - 5_000), n),
                ]
            )
        )
        p_plot = p[idx]
        exp = -np.log10((idx + 1) / (n + 1))
    else:
        p_plot = p
        exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs = -np.log10(p_plot)
    ax.scatter(exp, obs, s=4, alpha=0.45, c="#1f4e79", rasterized=True)
    m = max(float(np.nanmax(exp)), float(np.nanmax(obs)))
    ax.plot([0, m], [0, m], ls="--", c="gray", lw=1)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(f"{title}\nλGC≈{lam:.3f}  n={n:,}")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return lam


def manhattan_genome(df: pd.DataFrame, pcol: str, title: str, out_png: Path) -> None:
    d = df.dropna(subset=[pcol]).copy()
    d = d[np.isfinite(d[pcol]) & (d[pcol] > 0)]
    fig, ax = plt.subplots(figsize=(12, 3.5))
    if d.empty:
        ax.set_title(title + " (empty)")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    chrom_col = "CHR" if "CHR" in d.columns else ("CHROM" if "CHROM" in d.columns else None)
    d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
    d["nlp"] = -np.log10(d[pcol].astype(float))
    d = d.dropna(subset=["POS", "nlp"])

    if chrom_col is None:
        ax.scatter(d["POS"], d["nlp"], s=3, alpha=0.45, c="#1f4e79", rasterized=True)
        ax.set_xlabel("Position")
    else:
        d["_chrom"] = d[chrom_col].astype(str)
        chroms = sorted(d["_chrom"].unique(), key=chrom_sort_key)
        offset = 0.0
        ticks = []
        colors = ["#1f4e79", "#6baed6"]
        for i, chrom in enumerate(chroms):
            sub = d.loc[d["_chrom"] == chrom]
            x = sub["POS"].to_numpy(dtype=float) + offset
            ax.scatter(
                x,
                sub["nlp"],
                s=3,
                alpha=0.45,
                c=colors[i % 2],
                rasterized=True,
                label=chrom if len(chroms) <= 8 else None,
            )
            mid = offset + 0.5 * (sub["POS"].max() - sub["POS"].min() + 1)
            ticks.append((mid, chrom.replace("chr", "")))
            offset += float(sub["POS"].max()) + 1e6
        ax.set_xticks([t[0] for t in ticks])
        ax.set_xticklabels([t[1] for t in ticks], rotation=90, fontsize=7)
        ax.set_xlabel("Chromosome")

    ax.axhline(-np.log10(5e-8), color="crimson", ls="--", lw=0.8, label="5e-8")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def phenotype_counts(pheno_cov: Path | None, phenotype: str) -> dict:
    out = {
        "n_samples": np.nan,
        "n_cases": np.nan,
        "n_controls": np.nan,
        "n_missing_pheno": np.nan,
    }
    if pheno_cov is None or not pheno_cov.exists():
        return out
    df = pd.read_csv(pheno_cov, sep="\t", usecols=lambda c: c in {"ID", phenotype})
    if phenotype not in df.columns:
        return out
    y = pd.to_numeric(df[phenotype], errors="coerce")
    out["n_missing_pheno"] = int(y.isna().sum())
    y = y.dropna()
    out["n_samples"] = int(len(y))
    out["n_cases"] = int((y == 1).sum())
    out["n_controls"] = int((y == 0).sum())
    return out


def summarize_one(
    path: Path,
    phenotype: str,
    out_dir: Path,
    named_dir: Path,
    pheno_cov: Path | None,
    top_n: int,
    p_threshold: float,
    pcol: str = "P",
    named_suffix: str = ".tractor_mix.tsv",
    extra_p_columns: list[str] | None = None,
) -> dict:
    safe = sanitize_pheno(phenotype)
    ph_dir = out_dir / safe
    ph_dir.mkdir(parents=True, exist_ok=True)

    # Phenotype-named association table for Terra / notebooks
    named_path = named_dir / f"{safe}{named_suffix}"
    shutil.copy2(path, named_path)

    df = pd.read_csv(path, sep="\t")
    ac_pass = ac_pass_mask(df, pcol=pcol)
    n_var = len(df)
    n_pass = int(ac_pass.sum())
    n_joint = int(df[pcol].notna().sum()) if pcol in df.columns else 0

    lam_all = float("nan")
    lam_pass = float("nan")
    n_gw_sig = 0
    if pcol in df.columns:
        lam_all = qq_plot(
            df[pcol].to_numpy(dtype=float),
            f"{safe} {pcol} (all)",
            ph_dir / "qq_joint_all.png",
        )
        lam_pass = qq_plot(
            df.loc[ac_pass, pcol].to_numpy(dtype=float),
            f"{safe} {pcol} (AC-pass)",
            ph_dir / "qq_joint_acpass.png",
        )
        manhattan_genome(
            df.loc[ac_pass],
            pcol,
            f"{safe} {pcol} (AC-pass)",
            ph_dir / "manhattan_joint.png",
        )
        n_gw_sig = int(
            (
                ac_pass
                & df[pcol].notna()
                & (pd.to_numeric(df[pcol], errors="coerce") < p_threshold)
            ).sum()
        )

    extras = extra_p_columns
    if extras is None:
        extras = default_extra_p_columns(df, pcol)
    extra_lams: dict[str, float] = {}
    for c in extras:
        if c not in df.columns:
            continue
        extra_lams[f"lambda_gc_{c}"] = qq_plot(
            df.loc[ac_pass, c].to_numpy(dtype=float),
            f"{safe} {c} (AC-pass)",
            ph_dir / f"qq_{c}.png",
        )

    for c in [c for c in df.columns if c.startswith("Pval_anc")]:
        qq_plot(
            df.loc[ac_pass, c].to_numpy(dtype=float),
            f"{safe} {c} (AC-pass)",
            ph_dir / f"qq_{c}.png",
        )

    if pcol in df.columns:
        top = (
            df.loc[ac_pass]
            .dropna(subset=[pcol])
            .assign(_p=lambda x: pd.to_numeric(x[pcol], errors="coerce"))
            .sort_values("_p")
            .drop(columns=["_p"])
            .head(top_n)
        )
    else:
        top = df.head(0)
    top_path = ph_dir / "top_hits.tsv"
    top.to_csv(top_path, sep="\t", index=False)

    counts = phenotype_counts(pheno_cov, phenotype)
    rec = {
        "phenotype": phenotype,
        "phenotype_safe": safe,
        "p_column": pcol,
        "result_file": named_path.name,
        "result_path": str(named_path),
        "n_variants": n_var,
        "n_ac_pass": n_pass,
        "n_with_joint_P": n_joint,
        "n_genomewide_sig": n_gw_sig,
        "p_threshold": p_threshold,
        "lambda_gc_joint_all": lam_all,
        "lambda_gc_joint_acpass": lam_pass,
        "n_samples": counts["n_samples"],
        "n_cases": counts["n_cases"],
        "n_controls": counts["n_controls"],
        "n_missing_pheno": counts["n_missing_pheno"],
        "qq_acpass": str(ph_dir / "qq_joint_acpass.png"),
        "manhattan": str(ph_dir / "manhattan_joint.png"),
        "top_hits": str(top_path),
    }
    rec.update(extra_lams)
    return rec


def write_phewas_hits(summary: pd.DataFrame, out_dir: Path, named_dir: Path) -> Path | None:
    """Stack genome-wide significant rows across phenotypes for a PheWAS table."""
    rows = []
    for _, rec in summary.iterrows():
        path = named_dir / rec["result_file"]
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        pcol = str(rec["p_column"]) if "p_column" in rec.index and pd.notna(rec["p_column"]) else "P"
        if pcol not in df.columns:
            continue
        ac_pass = ac_pass_mask(df, pcol=pcol)
        thr = float(rec["p_threshold"])
        hit = df.loc[
            ac_pass
            & df[pcol].notna()
            & (pd.to_numeric(df[pcol], errors="coerce") < thr)
        ].copy()
        if hit.empty:
            continue
        hit.insert(0, "phenotype", rec["phenotype"])
        rows.append(hit)
    out = out_dir / "phewas_genomewide_hits.tsv"
    if not rows:
        pd.DataFrame({"phenotype": []}).to_csv(out, sep="\t", index=False)
        return out
    pd.concat(rows, ignore_index=True).to_csv(out, sep="\t", index=False)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--results",
        nargs="+",
        type=Path,
        required=True,
        help="Merged Tractor-Mix TSVs (parallel to --phenotypes if given).",
    )
    p.add_argument(
        "--phenotypes",
        nargs="*",
        default=None,
        help="Phenotype names parallel to --results (default: infer from filenames).",
    )
    p.add_argument("--pheno-cov", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--p-threshold", type=float, default=5e-8)
    p.add_argument(
        "--p-column",
        default="P",
        help="Primary p-value column (Tractor-Mix: P; FELIX: P_cct_admixed_c).",
    )
    p.add_argument(
        "--named-suffix",
        default=".tractor_mix.tsv",
        help="Suffix for phenotype-named copies under results_by_phenotype/.",
    )
    p.add_argument(
        "--extra-p-columns",
        nargs="*",
        default=None,
        help="Additional p-value columns for QQ plots. Default: FELIX _c_ extras "
        "when --p-column is P_cct_admixed_c; otherwise none besides Pval_anc*.",
    )
    p.add_argument(
        "--report-title",
        default="Tractor-Mix genome QC summary",
        help="Markdown heading for calibration_summary.md.",
    )
    args = p.parse_args()

    results = list(args.results)
    if args.phenotypes:
        if len(args.phenotypes) != len(results):
            raise SystemExit(
                f"--phenotypes length {len(args.phenotypes)} != --results {len(results)}"
            )
        phenos = list(args.phenotypes)
    else:
        phenos = []
        extra_suf = args.named_suffix if args.named_suffix.startswith(".") else f".{args.named_suffix}"
        for path in results:
            name = path.name
            for suf in (extra_suf, ".tractor_mix.tsv", ".felix.tsv", ".tsv"):
                if name.endswith(suf):
                    name = name[: -len(suf)]
                    break
            phenos.append(name)

    out_dir = args.out_dir
    named_dir = out_dir / "results_by_phenotype"
    qc_dir = out_dir / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    named_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path, pheno in zip(results, phenos):
        print(f"Summarizing {pheno} <- {path}")
        rows.append(
            summarize_one(
                path=path,
                phenotype=pheno,
                out_dir=qc_dir,
                named_dir=named_dir,
                pheno_cov=args.pheno_cov,
                top_n=args.top_n,
                p_threshold=args.p_threshold,
                pcol=args.p_column,
                named_suffix=args.named_suffix,
                extra_p_columns=args.extra_p_columns,
            )
        )

    summary = pd.DataFrame(rows)
    summary_path = out_dir / "calibration_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    manifest = summary[
        [
            "phenotype",
            "phenotype_safe",
            "result_file",
            "n_variants",
            "n_ac_pass",
            "lambda_gc_joint_acpass",
            "n_cases",
            "n_controls",
            "n_genomewide_sig",
        ]
    ].copy()
    manifest_path = out_dir / "results_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    # Wide λGC convenience table (single model run → one column)
    lam = summary.set_index("phenotype")[["lambda_gc_joint_acpass"]].rename(
        columns={"lambda_gc_joint_acpass": "lambda_gc"}
    )
    lam.to_csv(out_dir / "lambda_gc_wide.tsv", sep="\t")

    write_phewas_hits(summary, out_dir, named_dir)

    try:
        body = summary[
            [
                "phenotype",
                "n_cases",
                "n_controls",
                "n_ac_pass",
                "lambda_gc_joint_acpass",
                "n_genomewide_sig",
            ]
        ].to_markdown(index=False)
    except Exception:
        body = "```\n" + summary.to_string(index=False) + "\n```"

    (out_dir / "calibration_summary.md").write_text(
        f"# {args.report_title}\n\n"
        f"Primary p-value column: `{args.p_column}`.\n"
        "Per-phenotype association tables: `results_by_phenotype/`.\n"
        "Figures and top hits: `qc/<phenotype>/`.\n\n"
        + body
        + "\n"
    )
    print(f"Wrote {summary_path}")
    print(f"Wrote {manifest_path}")
    print(f"Wrote phenotype-named TSVs under {named_dir}")


if __name__ == "__main__":
    main()
