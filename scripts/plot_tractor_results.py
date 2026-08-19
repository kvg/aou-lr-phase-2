#!/usr/bin/env python3
"""QQ / Manhattan / AC-filter summaries for Tractor-Mix pilot result TSVs."""

from __future__ import annotations

import argparse
import math
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
    """Genomic control lambda using 1-df chi-square transform of p-values."""
    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size < 10:
        return float("nan")
    x = 2.0 * (erfcinv(p) ** 2)
    return float(np.median(x) / 0.4549364)


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
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs = -np.log10(p)
    ax.scatter(exp, obs, s=6, alpha=0.5, c="#1f4e79")
    m = max(float(exp.max()), float(obs.max()))
    ax.plot([0, m], [0, m], ls="--", c="gray", lw=1)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(f"{title}\nλGC≈{lam:.3f}  n={n:,}")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return lam


def manhattan(df: pd.DataFrame, pcol: str, title: str, out_png: Path) -> None:
    d = df.dropna(subset=[pcol]).copy()
    d = d[np.isfinite(d[pcol]) & (d[pcol] > 0)]
    fig, ax = plt.subplots(figsize=(10, 3.5))
    if d.empty:
        ax.set_title(title + " (empty)")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return
    d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
    d["nlp"] = -np.log10(d[pcol].astype(float))
    d = d.dropna(subset=["POS", "nlp"])
    ax.scatter(d["POS"], d["nlp"], s=4, alpha=0.5, c="#1f4e79")
    ax.axhline(-np.log10(5e-8), color="crimson", ls="--", lw=0.8, label="5e-8")
    ax.set_xlabel("Position")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def summarize_one(path: Path, out_dir: Path) -> dict:
    name = path.name
    pheno = name[: -len(".tsv")] if name.endswith(".tsv") else name
    df = pd.read_csv(path, sep="\t")
    out_ph = out_dir / pheno
    out_ph.mkdir(parents=True, exist_ok=True)

    include_cols = [c for c in df.columns if c.startswith("include_anc")]
    if include_cols:
        inc = df[include_cols].astype(str).apply(
            lambda s: s.str.lower().isin(["true", "1", "1.0", "t"])
        )
        ac_pass = inc.any(axis=1)
    else:
        ac_pass = df["P"].notna() if "P" in df.columns else pd.Series(False, index=df.index)

    n_var = len(df)
    n_pass = int(ac_pass.sum())
    n_joint = int(df["P"].notna().sum()) if "P" in df.columns else 0

    lam_all = float("nan")
    lam_pass = float("nan")
    if "P" in df.columns:
        lam_all = qq_plot(
            df["P"].to_numpy(dtype=float),
            f"{pheno} joint P (all rows with P)",
            out_ph / "qq_joint_all.png",
        )
        lam_pass = qq_plot(
            df.loc[ac_pass, "P"].to_numpy(dtype=float),
            f"{pheno} joint P (AC-pass)",
            out_ph / "qq_joint_acpass.png",
        )
        manhattan(
            df.loc[ac_pass], "P", f"{pheno} joint (AC-pass)", out_ph / "manhattan_joint.png"
        )

    for c in [c for c in df.columns if c.startswith("Pval_anc")]:
        qq_plot(df[c].to_numpy(dtype=float), f"{pheno} {c}", out_ph / f"qq_{c}.png")

    top = (
        df.dropna(subset=["P"]).sort_values("P").head(50) if "P" in df.columns else df.head(0)
    )
    top_path = out_ph / "top_hits.tsv"
    top.to_csv(top_path, sep="\t", index=False)

    return {
        "phenotype_file": str(path),
        "phenotype": pheno,
        "n_variants": n_var,
        "n_ac_pass": n_pass,
        "n_with_joint_P": n_joint,
        "lambda_gc_joint_all": lam_all,
        "lambda_gc_joint_acpass": lam_pass,
        "top_hits": str(top_path),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", nargs="+", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in args.results:
        print(f"Summarizing {path}")
        rows.append(summarize_one(path, args.out_dir))

    summary = pd.DataFrame(rows)
    summary_path = args.out_dir / "pilot_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    try:
        body = summary.to_markdown(index=False)
    except Exception:
        body = "```\n" + summary.to_string(index=False) + "\n```"
    (args.out_dir / "pilot_summary.md").write_text(
        "# Tractor-Mix pilot QC summary\n\n" + body + "\n"
    )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
