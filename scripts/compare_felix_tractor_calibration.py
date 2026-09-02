#!/usr/bin/env python3
"""Chr22 limited/full 2×2 calibration gate: FELIX vs Tractor-Mix.

Expects up to four labeled result directories sharing phenotypes:
  felix_limited, felix_full, tractor_limited, tractor_full

Writes combined QQ plots, λGC summary (calibration_summary.tsv), and flags
sample-count mismatches as comparability failures.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MODEL_ORDER = ["felix_limited", "felix_full", "tractor_limited", "tractor_full"]
MODEL_COLORS = {
    "felix_limited": "#9467bd",
    "felix_full": "#8c564b",
    "tractor_limited": "#1f77b4",
    "tractor_full": "#ff7f0e",
}

FELIX_P_COLS = [
    "P_cct_admixed_c",
    "P_cct_admixed",
    "P_hom_admixed_c",
    "P_het_admixed_c",
    "p.value_ancALL",
]
TRACTOR_P_COLS = ["P", "pvalue", "p.value"]


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


def phenotype_key(path: Path) -> str:
    name = path.name
    for suffix in [
        ".felix.tsv",
        ".tractor_mix.tsv",
        ".tractor_mix.txt",
        ".tsv",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def discover_results(directory: Path) -> dict[str, Path]:
    if not directory.exists():
        return {}
    out: dict[str, Path] = {}
    for path in sorted(directory.glob("*.tsv")):
        if path.name.endswith(".raw.txt") or "summary" in path.name:
            continue
        out[phenotype_key(path)] = path
    return out


def pick_p_column(df: pd.DataFrame, model: str) -> str | None:
    cols = FELIX_P_COLS if model.startswith("felix") else FELIX_P_COLS + TRACTOR_P_COLS
    for c in cols:
        if c in df.columns:
            return c
    return None


def load_pvalues(path: Path, model: str) -> tuple[pd.DataFrame, np.ndarray, dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    meta: dict = {"n_rows": len(df)}
    pcol = pick_p_column(df, model)
    if pcol is None:
        raise SystemExit(f"No p-value column in {path}; columns={list(df.columns)[:30]}")

    n_eff = None
    for c in ["n_samples", "N", "N_eff", "Neff"]:
        if c in df.columns:
            vals = pd.to_numeric(df[c], errors="coerce").dropna()
            if not vals.empty:
                n_eff = float(vals.median())
                break
    meta["n_eff"] = n_eff

    include_cols = [c for c in df.columns if c.startswith("include_anc")]
    if include_cols:
        inc = df[include_cols].astype(str).apply(
            lambda s: s.str.lower().isin(["true", "1", "1.0", "t"])
        )
        keep = inc.any(axis=1)
        meta["filter"] = "tractor_ac_pass"
    else:
        keep = df[pcol].notna()
        meta["filter"] = "pvalue_present"

    sub = df.loc[keep].copy()
    p = pd.to_numeric(sub[pcol], errors="coerce").to_numpy(dtype=float)
    meta["n_tested"] = int(np.isfinite(p).sum())
    meta["pcol"] = pcol

    if {"CHR", "POS"}.issubset(sub.columns):
        sub["_varkey"] = sub["CHR"].astype(str) + ":" + sub["POS"].astype(str)
    elif "MarkerID" in sub.columns:
        sub["_varkey"] = sub["MarkerID"].astype(str)
    else:
        sub["_varkey"] = np.arange(len(sub)).astype(str)

    sub["_p"] = pd.to_numeric(sub[pcol], errors="coerce")
    return sub, p, meta


def combined_qq(series: dict[str, np.ndarray], title: str, out_png: Path) -> dict[str, float]:
    fig, ax = plt.subplots(figsize=(6, 6))
    lams: dict[str, float] = {}
    max_m = 0.0
    for model in MODEL_ORDER:
        if model not in series:
            continue
        p = np.sort(series[model][np.isfinite(series[model]) & (series[model] > 0) & (series[model] <= 1)])
        lam = lambda_gc_from_p(p)
        lams[model] = lam
        if p.size == 0:
            continue
        n = p.size
        exp = -np.log10(np.arange(1, n + 1) / (n + 1))
        obs = -np.log10(p)
        max_m = max(max_m, float(exp.max()), float(obs.max()))
        ax.scatter(
            exp,
            obs,
            s=5,
            alpha=0.45,
            c=MODEL_COLORS[model],
            label=f"{model} λ={lam:.3f} n={n:,}",
        )
    if max_m > 0:
        ax.plot([0, max_m], [0, max_m], ls="--", c="gray", lw=1)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return lams


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--felix-limited-dir", type=Path, default=None)
    p.add_argument("--felix-full-dir", type=Path, default=None)
    p.add_argument("--tractor-limited-dir", type=Path, default=None)
    p.add_argument("--tractor-full-dir", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    model_dirs = {
        "felix_limited": args.felix_limited_dir,
        "felix_full": args.felix_full_dir,
        "tractor_limited": args.tractor_limited_dir,
        "tractor_full": args.tractor_full_dir,
    }
    model_dirs = {k: v for k, v in model_dirs.items() if v is not None}

    discovered = {m: discover_results(d) for m, d in model_dirs.items()}
    phenotypes = sorted({ph for d in discovered.values() for ph in d})
    if not phenotypes:
        raise SystemExit("No phenotype result TSVs found in provided directories")

    summary_rows = []
    comparability_flags = []

    for pheno in phenotypes:
        ph_dir = args.out_dir / pheno
        ph_dir.mkdir(parents=True, exist_ok=True)
        p_series: dict[str, np.ndarray] = {}
        n_effs: dict[str, float | None] = {}

        for model, files in discovered.items():
            if pheno not in files:
                continue
            _df, pvals, meta = load_pvalues(files[pheno], model)
            p_series[model] = pvals
            n_effs[model] = meta.get("n_eff")
            summary_rows.append(
                {
                    "phenotype": pheno,
                    "model": model,
                    "result_file": str(files[pheno]),
                    "p_column": meta["pcol"],
                    "n_variants_tested": meta["n_tested"],
                    "n_eff": meta.get("n_eff"),
                    "lambda_gc": lambda_gc_from_p(pvals),
                    "filter": meta["filter"],
                }
            )

        combined_qq(p_series, f"{pheno} FELIX vs Tractor λGC", ph_dir / "qq_matched.png")

        present_ns = {m: n for m, n in n_effs.items() if n is not None}
        if len(present_ns) >= 2:
            vals = list(present_ns.values())
            if max(vals) != min(vals):
                msg = (
                    f"COMPARABILITY FAILURE for {pheno}: sample counts differ "
                    f"across models: {present_ns}"
                )
                comparability_flags.append(
                    {"phenotype": pheno, "status": "FAIL", "detail": msg}
                )
                (ph_dir / "COMPARABILITY_FAIL.txt").write_text(msg + "\n")
            else:
                comparability_flags.append(
                    {"phenotype": pheno, "status": "OK", "detail": f"n_eff={vals[0]}"}
                )

    summary = pd.DataFrame(summary_rows)
    summary_path = args.out_dir / "calibration_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    if not summary.empty:
        wide = summary.pivot_table(
            index="phenotype", columns="model", values="lambda_gc", aggfunc="first"
        )
        wide.to_csv(args.out_dir / "lambda_gc_wide.tsv", sep="\t")

    pd.DataFrame(comparability_flags).to_csv(
        args.out_dir / "comparability_flags.tsv", sep="\t", index=False
    )

    lines = [
        "# FELIX vs Tractor-Mix chr22 calibration gate",
        "",
        "Four-way comparison: limited/full covariates × FELIX / Tractor-Mix.",
        "Target: λGC ≈ 1 on null phenotypes; treat sample-count mismatches as FAIL.",
        "",
    ]
    if not summary.empty:
        try:
            lines.append(summary.to_markdown(index=False))
        except Exception:
            lines.append("```\n" + summary.to_string(index=False) + "\n```")
    (args.out_dir / "calibration_summary.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {summary_path} ({len(phenotypes)} phenotypes, models={list(model_dirs)})")


if __name__ == "__main__":
    main()
