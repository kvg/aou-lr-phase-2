#!/usr/bin/env python3
"""Matched four-model calibration QC for Tractor-Mix and SAIGE pilots.

Expects up to four labeled result directories (or file globs) sharing phenotypes:
  tractor_limited, tractor_full, saige_limited, saige_full

Writes combined QQ plots, λGC / N summary, and secondary top-hit overlap.
Flags sample-count mismatches as comparability failures rather than
interpreting λGC differences.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MODEL_ORDER = ["tractor_limited", "tractor_full", "saige_limited", "saige_full"]
MODEL_COLORS = {
    "tractor_limited": "#1f77b4",
    "tractor_full": "#ff7f0e",
    "saige_limited": "#2ca02c",
    "saige_full": "#d62728",
}


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
        ".tractor_mix.tsv",
        ".saige.tsv",
        ".tractor_mix.txt",
        ".saige.txt",
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


def load_pvalues(path: Path) -> tuple[pd.DataFrame, np.ndarray, dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    meta: dict = {"n_rows": len(df)}

    # Tractor-Mix joint P vs SAIGE normalized pvalue
    pcol = None
    for c in ["P", "pvalue", "p.value", "Pvalue"]:
        if c in df.columns:
            pcol = c
            break
    if pcol is None:
        raise SystemExit(f"No p-value column in {path}; columns={list(df.columns)[:30]}")

    # Effective sample size if present in companion meta or columns
    n_eff = None
    for c in ["n_samples", "N", "N_eff", "Neff"]:
        if c in df.columns:
            vals = pd.to_numeric(df[c], errors="coerce").dropna()
            if not vals.empty:
                n_eff = float(vals.median())
                break
    meta["n_eff"] = n_eff

    # AC-pass filter for Tractor; MAC filter already applied in SAIGE step2
    include_cols = [c for c in df.columns if c.startswith("include_anc")]
    if include_cols:
        inc = df[include_cols].astype(str).apply(
            lambda s: s.str.lower().isin(["true", "1", "1.0", "t"])
        )
        keep = inc.any(axis=1)
        meta["filter"] = "tractor_ac_pass"
    elif "mac" in df.columns:
        mac = pd.to_numeric(df["mac"], errors="coerce")
        keep = mac.notna()
        meta["filter"] = "saige_mac_present"
    else:
        keep = df[pcol].notna()
        meta["filter"] = "pvalue_present"

    sub = df.loc[keep].copy()
    p = pd.to_numeric(sub[pcol], errors="coerce").to_numpy(dtype=float)
    meta["n_tested"] = int(np.isfinite(p).sum())
    meta["pcol"] = pcol

    # Variant key for overlap
    if {"CHR", "POS"}.issubset(sub.columns) or {"chrom", "pos"}.issubset(sub.columns):
        chrom = sub["CHR"] if "CHR" in sub.columns else sub["chrom"]
        pos = sub["POS"] if "POS" in sub.columns else sub["pos"]
        sub["_varkey"] = chrom.astype(str) + ":" + pos.astype(str)
    elif "marker_id" in sub.columns:
        sub["_varkey"] = sub["marker_id"].astype(str)
    elif "SNPID" in sub.columns:
        sub["_varkey"] = sub["SNPID"].astype(str)
    else:
        sub["_varkey"] = np.arange(len(sub)).astype(str)

    sub["_p"] = pd.to_numeric(sub[pcol], errors="coerce")
    return sub, p, meta


def load_null_meta(meta_dir: Path | None, phenotype: str) -> int | None:
    if meta_dir is None or not meta_dir.exists():
        return None
    candidates = list(meta_dir.glob(f"{phenotype}*.null_meta.tsv")) + list(
        meta_dir.glob(f"*{phenotype}*.tsv")
    )
    for path in candidates:
        try:
            m = pd.read_csv(path, sep="\t")
            if "n_samples" in m.columns:
                return int(m["n_samples"].iloc[0])
        except Exception:
            continue
    return None


def combined_qq(
    series: dict[str, np.ndarray],
    title: str,
    out_png: Path,
) -> dict[str, float]:
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


def top_overlap(frames: dict[str, pd.DataFrame], k: int = 50) -> pd.DataFrame:
    tops: dict[str, set[str]] = {}
    for model, df in frames.items():
        d = df.dropna(subset=["_p"]).sort_values("_p").head(k)
        tops[model] = set(d["_varkey"].astype(str))
    rows = []
    models = [m for m in MODEL_ORDER if m in tops]
    for i, a in enumerate(models):
        for b in models[i:]:
            inter = len(tops[a] & tops[b])
            rows.append(
                {
                    "model_a": a,
                    "model_b": b,
                    "k": k,
                    "overlap": inter,
                    "jaccard": inter / max(1, len(tops[a] | tops[b])),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tractor-limited-dir", type=Path, default=None)
    p.add_argument("--tractor-full-dir", type=Path, default=None)
    p.add_argument("--saige-limited-dir", type=Path, default=None)
    p.add_argument("--saige-full-dir", type=Path, default=None)
    p.add_argument(
        "--null-meta-dirs",
        nargs="*",
        default=[],
        help="Optional dirs with *.null_meta.tsv for effective N (SAIGE).",
    )
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--top-k", type=int, default=50)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    model_dirs = {
        "tractor_limited": args.tractor_limited_dir,
        "tractor_full": args.tractor_full_dir,
        "saige_limited": args.saige_limited_dir,
        "saige_full": args.saige_full_dir,
    }
    model_dirs = {k: v for k, v in model_dirs.items() if v is not None}

    discovered = {m: discover_results(d) for m, d in model_dirs.items()}
    phenotypes = sorted({ph for d in discovered.values() for ph in d})
    if not phenotypes:
        raise SystemExit("No phenotype result TSVs found in provided directories")

    null_meta_dirs = [Path(x) for x in args.null_meta_dirs]

    summary_rows = []
    overlap_rows = []
    comparability_flags = []

    for pheno in phenotypes:
        ph_dir = args.out_dir / pheno
        ph_dir.mkdir(parents=True, exist_ok=True)
        p_series: dict[str, np.ndarray] = {}
        frames: dict[str, pd.DataFrame] = {}
        n_effs: dict[str, float | None] = {}

        for model, files in discovered.items():
            if pheno not in files:
                continue
            df, pvals, meta = load_pvalues(files[pheno])
            frames[model] = df
            p_series[model] = pvals
            n_from_meta = None
            for md in null_meta_dirs:
                n_from_meta = load_null_meta(md, pheno)
                if n_from_meta is not None:
                    break
            n_eff = meta.get("n_eff") if meta.get("n_eff") is not None else n_from_meta
            n_effs[model] = n_eff
            summary_rows.append(
                {
                    "phenotype": pheno,
                    "model": model,
                    "result_file": str(files[pheno]),
                    "n_variants_tested": meta["n_tested"],
                    "n_eff": n_eff,
                    "lambda_gc": lambda_gc_from_p(pvals),
                    "filter": meta["filter"],
                }
            )

        lams = combined_qq(p_series, f"{pheno} matched QQ", ph_dir / "qq_matched.png")
        for model, lam in lams.items():
            # already in summary_rows; ensure consistency
            pass

        # Comparability: sample counts must match across models when available
        present_ns = {m: n for m, n in n_effs.items() if n is not None}
        if len(present_ns) >= 2:
            vals = list(present_ns.values())
            if max(vals) != min(vals):
                msg = (
                    f"COMPARABILITY FAILURE for {pheno}: sample counts differ "
                    f"across models: {present_ns}. Do not interpret λGC differences."
                )
                comparability_flags.append(
                    {"phenotype": pheno, "status": "FAIL", "detail": msg, **{f"n_{k}": v for k, v in present_ns.items()}}
                )
                (ph_dir / "COMPARABILITY_FAIL.txt").write_text(msg + "\n")
            else:
                comparability_flags.append(
                    {
                        "phenotype": pheno,
                        "status": "OK",
                        "detail": f"n_eff matched at {vals[0]}",
                        **{f"n_{k}": v for k, v in present_ns.items()},
                    }
                )
        else:
            comparability_flags.append(
                {
                    "phenotype": pheno,
                    "status": "UNKNOWN",
                    "detail": "n_eff not available for >=2 models; rely on shared analysis_samples prep",
                }
            )

        if len(frames) >= 2:
            ov = top_overlap(frames, k=args.top_k)
            ov.insert(0, "phenotype", pheno)
            overlap_rows.append(ov)
            ov.to_csv(ph_dir / "top_hit_overlap.tsv", sep="\t", index=False)

            # Secondary top-hit table union
            parts = []
            for model, df in frames.items():
                top = df.dropna(subset=["_p"]).sort_values("_p").head(args.top_k).copy()
                top.insert(0, "model", model)
                parts.append(top)
            pd.concat(parts, ignore_index=True).to_csv(
                ph_dir / "top_hits_by_model.tsv", sep="\t", index=False
            )

    summary = pd.DataFrame(summary_rows)
    summary_path = args.out_dir / "calibration_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    # Wide λGC table
    if not summary.empty:
        wide = summary.pivot_table(
            index="phenotype", columns="model", values="lambda_gc", aggfunc="first"
        )
        wide.to_csv(args.out_dir / "lambda_gc_wide.tsv", sep="\t")
        nwide = summary.pivot_table(
            index="phenotype", columns="model", values="n_variants_tested", aggfunc="first"
        )
        nwide.to_csv(args.out_dir / "n_tested_wide.tsv", sep="\t")

    flags = pd.DataFrame(comparability_flags)
    flags.to_csv(args.out_dir / "comparability_flags.tsv", sep="\t", index=False)

    if overlap_rows:
        pd.concat(overlap_rows, ignore_index=True).to_csv(
            args.out_dir / "top_hit_overlap_all.tsv", sep="\t", index=False
        )

    # Markdown report
    lines = [
        "# Matched four-model calibration QC",
        "",
        "Models: Tractor-Mix / SAIGE × limited / full covariates on the shared cohort.",
        "",
        "## Comparability",
        "",
        "Sample-count differences are treated as **comparability failures**; "
        "do not interpret λGC differences when status=FAIL.",
        "",
    ]
    if not flags.empty:
        try:
            lines.append(flags.to_markdown(index=False))
        except Exception:
            lines.append("```\n" + flags.to_string(index=False) + "\n```")
    lines.extend(["", "## λGC and tested-variant counts", ""])
    if not summary.empty:
        try:
            lines.append(summary.to_markdown(index=False))
        except Exception:
            lines.append("```\n" + summary.to_string(index=False) + "\n```")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Combined QQ plots: `<phenotype>/qq_matched.png`",
            "- Top-hit overlap is secondary to calibration.",
            "",
        ]
    )
    (args.out_dir / "calibration_summary.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {summary_path}")
    print(f"Phenotypes: {len(phenotypes)}; models: {list(model_dirs)}")


if __name__ == "__main__":
    main()
