#!/usr/bin/env python3
"""Tractor null-λ pilot for FLARE recipe survivors (decisive metric, Part 4).

1. Slice recipe ``anc.vcf`` to the evaluation window.
2. ``extract-tracts-flare --num-ancs 5`` on analysis samples.
3. For each of K phenotype permutations: shuffle a binary phenotype column,
   ``fit_null.R`` (reuse GRM), ``tractor-mix-score`` with lowered AC threshold,
   compute λGC.
4. Aggregate mean |λ−1| and bootstrap CI; declare winner among survivors.
   Tie-break overlapping CIs with Part 2 ``mean_ll`` only within that set.

If fewer than ``min_tested_sites`` (default 100) remain after the AC filter,
widen/lower AC and record ``unstable_lambda=true`` — do not silently compare.

Example::

    python3 scripts/flare_lai_null_lambda.py \\
      --anc-vcf recipe.anc.vcf.gz \\
      --region chr22:26897597-36897597 \\
      --samples analysis_samples.txt \\
      --pheno-cov pheno_cov.tsv \\
      --phenotype some_binary \\
      --covariates age,sex,PC1,PC2,PC3,PC4,PC5 \\
      --grm-rds sparse_grm.rds \\
      --fit-null-r scripts/fit_null.R \\
      --out-dir null_lambda/recipe_x \\
      --n-perm 50 --seed 1 --ac-threshold 20
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


def _pd():
    import pandas as pd

    return pd


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
    """Genomic control λ (1-df χ²); same formula as plot_tractor_results.py."""
    p = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    if p.size < 10:
        return float("nan")
    x = 2.0 * (erfcinv(p) ** 2)
    return float(np.median(x) / 0.4549364)


def bootstrap_mean_ci(
    values: Sequence[float],
    *,
    n_boot: int = 1000,
    seed: int = 1,
    alpha: float = 0.05,
) -> dict[str, float]:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return {"mean": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        means[i] = float(rng.choice(arr, size=arr.size, replace=True).mean())
    lo = float(np.quantile(means, alpha / 2))
    hi = float(np.quantile(means, 1 - alpha / 2))
    return {"mean": float(arr.mean()), "ci_low": lo, "ci_high": hi}


def abs_dev_from_one(lambdas: Sequence[float]) -> list[float]:
    return [abs(float(x) - 1.0) for x in lambdas if np.isfinite(x)]


def cis_overlap(a: dict[str, float], b: dict[str, float]) -> bool:
    if not all(np.isfinite([a["ci_low"], a["ci_high"], b["ci_low"], b["ci_high"]])):
        return False
    return not (a["ci_high"] < b["ci_low"] or b["ci_high"] < a["ci_low"])


def pick_winner(
    rows: list[dict[str, Any]],
    *,
    lambda_key: str = "abs_lambda_dev",
    ll_key: str = "mean_ll",
) -> dict[str, Any]:
    """Winner = lowest mean |λ−1|; tie-break overlapping CIs with higher mean_ll."""
    usable = [r for r in rows if not r.get("unstable_lambda") and np.isfinite(r.get(lambda_key, float("nan")))]
    if not usable:
        return {"winner": None, "reason": "no_stable_survivors", "tied": []}
    usable.sort(key=lambda r: float(r[lambda_key]))
    best = usable[0]
    tied = [
        r
        for r in usable
        if cis_overlap(best["abs_lambda_ci"], r["abs_lambda_ci"])
    ]
    if len(tied) == 1:
        return {"winner": best.get("experiment"), "reason": "unique_min_abs_lambda", "tied": [best.get("experiment")]}
    # tie-break on mean_ll among overlapping set
    with_ll = [r for r in tied if np.isfinite(r.get(ll_key, float("nan")))]
    if not with_ll:
        return {
            "winner": best.get("experiment"),
            "reason": "ci_overlap_no_ll",
            "tied": [r.get("experiment") for r in tied],
        }
    with_ll.sort(key=lambda r: -float(r[ll_key]))
    return {
        "winner": with_ll[0].get("experiment"),
        "reason": "ci_overlap_mean_ll_tiebreak",
        "tied": [r.get("experiment") for r in tied],
    }


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, text=True, **kwargs)


def slice_anc_vcf(anc_vcf: str, region: str, out_vcf: Path) -> Path:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "bcftools",
            "view",
            "-r",
            region,
            "-Oz",
            "-o",
            str(out_vcf),
            anc_vcf,
        ]
    )
    run(["bcftools", "index", "-f", str(out_vcf)])
    return out_vcf


def extract_tracts(
    anc_vcf: Path,
    samples: Path,
    output_dir: Path,
    *,
    num_ancs: int = 5,
    extract_bin: str = "extract-tracts-flare",
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    run(
        [
            extract_bin,
            "--vcf",
            str(anc_vcf),
            "--num-ancs",
            str(num_ancs),
            "--output-dir",
            str(output_dir),
            "--compress-output",
            "--samples",
            str(samples),
        ]
    )
    # Order like TractorMixPilot: anc_00.dosage.txt.gz ...
    ordered_dir = output_dir / "ordered"
    ordered_dir.mkdir(parents=True, exist_ok=True)
    dosage_files: list[Path] = []
    for p in sorted(output_dir.glob("*.dosage.txt.gz")):
        if ".collapse." in p.name or ".split." in p.name:
            continue
        m = re.search(r"(?<![.])\.anc(\d+)\.dosage\.txt\.gz$", p.name)
        if not m:
            # also accept basename.ancN.dosage.txt.gz
            m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
            if not m or ".collapse." in p.name or ".split." in p.name:
                continue
        dest = ordered_dir / f"anc_{int(m.group(1)):02d}.dosage.txt.gz"
        shutil.copy2(p, dest)
        dosage_files.append(dest)
    dosage_files = sorted(set(dosage_files))
    if len(dosage_files) != num_ancs:
        raise SystemExit(
            f"expected {num_ancs} dosage files, found {len(dosage_files)} in {output_dir}"
        )
    return dosage_files


def shuffle_pheno(
    pheno_cov: Path,
    phenotype: str,
    out_path: Path,
    *,
    seed: int,
    id_col: str = "IID",
) -> Path:
    df = _pd().read_csv(pheno_cov, sep="\t", dtype=str)
    if phenotype not in df.columns:
        # try common alternatives
        raise SystemExit(f"phenotype {phenotype!r} not in {pheno_cov}")
    rng = random.Random(seed)
    vals = list(df[phenotype])
    rng.shuffle(vals)
    df[phenotype] = vals
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def fit_null(
    *,
    fit_null_r: Path,
    pheno_cov: Path,
    phenotype: str,
    covariates: str,
    grm_rds: Path,
    out_null_rds: Path,
    out_null_export: Path,
) -> None:
    run(
        [
            "Rscript",
            str(fit_null_r),
            "--pheno-cov",
            str(pheno_cov),
            "--phenotype",
            phenotype,
            "--covariates",
            covariates,
            "--grm-rds",
            str(grm_rds),
            "--out-null-rds",
            str(out_null_rds),
            "--out-null-export",
            str(out_null_export),
        ]
    )


def run_score(
    *,
    score_bin: str,
    null_export: Path,
    dosage_files: Sequence[Path],
    out_tsv: Path,
    ac_threshold: int,
    threads: int = 8,
) -> Path:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        score_bin,
        "--null-export",
        str(null_export),
        "--dosage-files",
        *[str(p) for p in dosage_files],
        "--ac-threshold",
        str(ac_threshold),
        "--threads",
        str(threads),
        "--out",
        str(out_tsv),
    ]
    run(cmd)
    return out_tsv


def read_pvalues(score_tsv: Path, p_col: str = "P") -> np.ndarray:
    df = _pd().read_csv(score_tsv, sep="\t")
    # tractor-mix-score column names vary; try common ones
    for cand in (p_col, "p", "Pvalue", "pvalue", "P_CCT", "p_CCT"):
        if cand in df.columns:
            return _pd().to_numeric(df[cand], errors="coerce").to_numpy(dtype=float)
    # last resort: any column starting with p
    for c in df.columns:
        if str(c).lower().startswith("p"):
            return _pd().to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
    raise SystemExit(f"no p-value column in {score_tsv}; columns={list(df.columns)}")


def aggregate_perm_lambdas(
    lambdas: Sequence[float],
    *,
    n_tested_sites: int,
    min_tested_sites: int,
    seed: int,
) -> dict[str, Any]:
    abs_devs = abs_dev_from_one(lambdas)
    ci = bootstrap_mean_ci(abs_devs, seed=seed)
    unstable = n_tested_sites < min_tested_sites or not np.isfinite(ci["mean"])
    return {
        "lambdas": [float(x) for x in lambdas],
        "n_perm": len(list(lambdas)),
        "mean_lambda": float(np.nanmean(lambdas)) if len(lambdas) else float("nan"),
        "abs_lambda_dev": ci["mean"],
        "abs_lambda_ci": ci,
        "n_tested_sites": n_tested_sites,
        "min_tested_sites": min_tested_sites,
        "unstable_lambda": bool(unstable),
    }


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--anc-vcf", required=True)
    p.add_argument("--region", required=True)
    p.add_argument("--samples", required=True, type=Path)
    p.add_argument("--pheno-cov", required=True, type=Path)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--covariates", required=True, help="comma-separated covariate names")
    p.add_argument("--grm-rds", required=True, type=Path)
    p.add_argument("--fit-null-r", type=Path, default=Path("scripts/fit_null.R"))
    p.add_argument("--extract-bin", default="extract-tracts-flare")
    p.add_argument("--score-bin", default="tractor-mix-score")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--n-perm", type=int, default=50)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ac-threshold", type=int, default=20)
    p.add_argument("--min-tested-sites", type=int, default=100)
    p.add_argument("--num-ancs", type=int, default=5)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--experiment", default="")
    p.add_argument(
        "--aggregate-only",
        type=Path,
        default=None,
        help="JSON list of λ values (+ n_tested_sites) to aggregate without running Tractor",
    )
    p.add_argument(
        "--declare-winner",
        type=Path,
        default=None,
        help="JSON list of per-recipe score dicts → pick winner",
    )
    p.add_argument("--dry-run", action="store_true", help="Write plan JSON only")
    args = p.parse_args(argv)

    if args.declare_winner is not None:
        rows = json.loads(args.declare_winner.read_text())
        decision = pick_winner(rows)
        out = args.out_dir / "winner.json"
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(decision, indent=2) + "\n")
        print(json.dumps(decision, indent=2))
        print("wrote", out)
        return 0

    if args.aggregate_only is not None:
        payload = json.loads(args.aggregate_only.read_text())
        lambdas = payload["lambdas"]
        n_sites = int(payload.get("n_tested_sites", 0))
        summary = aggregate_perm_lambdas(
            lambdas,
            n_tested_sites=n_sites,
            min_tested_sites=args.min_tested_sites,
            seed=args.seed,
        )
        summary["experiment"] = args.experiment or payload.get("experiment", "")
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out = args.out_dir / "null_lambda_score.json"
        out.write_text(json.dumps(summary, indent=2) + "\n")
        print("wrote", out)
        return 0

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    plan = {
        "experiment": args.experiment,
        "anc_vcf": args.anc_vcf,
        "region": args.region,
        "n_perm": args.n_perm,
        "seed": args.seed,
        "ac_threshold": args.ac_threshold,
        "min_tested_sites": args.min_tested_sites,
        "steps": [
            "bcftools view -r REGION anc.vcf",
            "extract-tracts-flare",
            "for k in 1..K: shuffle pheno → fit_null.R → tractor-mix-score → λGC",
            "bootstrap CI on |λ-1|; flag unstable if n_tested_sites < min",
        ],
    }
    (out_dir / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    if args.dry_run:
        print("wrote", out_dir / "plan.json")
        return 0

    for bin_name in (args.extract_bin, args.score_bin):
        if shutil.which(bin_name) is None:
            raise SystemExit(
                f"{bin_name} not found on PATH; use --dry-run or --aggregate-only for offline scoring"
            )

    sliced = slice_anc_vcf(args.anc_vcf, args.region, out_dir / "window.anc.vcf.gz")
    extract_dir = out_dir / "extract"
    dosage_files = extract_tracts(
        sliced,
        args.samples,
        extract_dir,
        num_ancs=args.num_ancs,
        extract_bin=args.extract_bin,
    )

    lambdas: list[float] = []
    n_tested_sites = 0
    for k in range(args.n_perm):
        seed_k = args.seed + k
        perm_dir = out_dir / f"perm_{k:03d}"
        perm_dir.mkdir(parents=True, exist_ok=True)
        pheno_k = shuffle_pheno(
            args.pheno_cov, args.phenotype, perm_dir / "pheno_cov.tsv", seed=seed_k
        )
        null_rds = perm_dir / "null.rds"
        null_export = perm_dir / "null_export"
        fit_null(
            fit_null_r=args.fit_null_r,
            pheno_cov=pheno_k,
            phenotype=args.phenotype,
            covariates=args.covariates,
            grm_rds=args.grm_rds,
            out_null_rds=null_rds,
            out_null_export=null_export,
        )
        score_tsv = perm_dir / "scores.tsv"
        run_score(
            score_bin=args.score_bin,
            null_export=null_export,
            dosage_files=dosage_files,
            out_tsv=score_tsv,
            ac_threshold=args.ac_threshold,
            threads=args.threads,
        )
        pvals = read_pvalues(score_tsv)
        n_tested_sites = max(n_tested_sites, int(np.isfinite(pvals).sum()))
        lam = lambda_gc_from_p(pvals)
        lambdas.append(lam)
        (perm_dir / "lambda.json").write_text(
            json.dumps({"perm": k, "lambda_gc": lam, "n_p": int(np.isfinite(pvals).sum())}, indent=2)
            + "\n"
        )

    summary = aggregate_perm_lambdas(
        lambdas,
        n_tested_sites=n_tested_sites,
        min_tested_sites=args.min_tested_sites,
        seed=args.seed,
    )
    summary.update(
        {
            "experiment": args.experiment,
            "region": args.region,
            "ac_threshold": args.ac_threshold,
            "phenotype": args.phenotype,
        }
    )
    if summary["unstable_lambda"]:
        summary["unstable_note"] = (
            f"n_tested_sites={n_tested_sites} < min_tested_sites={args.min_tested_sites}; "
            "widen window or lower --ac-threshold before comparing λ"
        )
    out = out_dir / "null_lambda_score.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
