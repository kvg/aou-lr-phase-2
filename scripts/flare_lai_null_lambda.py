#!/usr/bin/env python3
"""Tractor null-λ pilot for FLARE recipe survivors (decisive metric, Part 4).

Phenotypes are **directly simulated** (not shuffled):

    L_i = β·q_i + γ_GC + u_family + ε

with fixed global ancestry ``q`` (same file for every recipe), GC fixed
effects, and a shared family random effect for the Phase-2 pedigrees.
Liability is thresholded to match an anchor phenotype's case rate.

**Preflights (fail closed):**
1. ``q`` vs null-model PCs — if multivariate R² is too high, drop PCs from
   this pilot's null (production Tractor null need not match).
2. Anchor phenotype must have a non-trivial ancestry log-OR before locking
   β_mid; otherwise try fallbacks.

Run β at lo/mid/hi (mid from validated anchor; ±2× bracket). Trust a winner
only if the same recipe wins at all three magnitudes.

Example::

    python3 scripts/flare_lai_null_lambda.py preflight \\
      --pheno-cov pheno_cov.tsv \\
      --global-anc fixed.global.anc.gz \\
      --covariates age,sex,GC,PC1,PC2,PC3,PC4,PC5 \\
      --anchor-candidates trait_a,trait_b \\
      --out null_lambda_preflight.json

    python3 scripts/flare_lai_null_lambda.py run \\
      --anc-vcf recipe.anc.vcf.gz --region chr22:26897597-36897597 \\
      --samples analysis_samples.txt \\
      --pheno-cov pheno_cov.tsv --global-anc fixed.global.anc.gz \\
      --preflight null_lambda_preflight.json \\
      --grm-rds sparse_grm.rds --fit-null-r scripts/fit_null.R \\
      --out-dir null_lambda/recipe_x --n-rep 20 --seed 1
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np

ANC_NAMES = ("eas", "amr", "eur", "afr", "sas")
DEFAULT_R2_MAX = 0.90
DEFAULT_OR_ABS_MIN = 0.3
DEFAULT_OR_PMAX = 0.01
SIM_PHENO_COL = "sim_pheno"


def _pd():
    import pandas as pd

    return pd


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


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
    preflight: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Winner = lowest mean |λ−1|; tie-break overlapping CIs with higher mean_ll."""
    if preflight is not None:
        gate = assert_preflight_ok(preflight)
        if not gate["ok"]:
            return {
                "winner": None,
                "reason": gate["reason"],
                "tied": [],
                "ranking_unstable": True,
                "preflight": gate,
            }
    usable = [
        r for r in rows if not r.get("unstable_lambda") and np.isfinite(r.get(lambda_key, float("nan")))
    ]
    if not usable:
        return {"winner": None, "reason": "no_stable_survivors", "tied": [], "ranking_unstable": True}
    usable.sort(key=lambda r: float(r[lambda_key]))
    best = usable[0]
    tied = [r for r in usable if cis_overlap(best["abs_lambda_ci"], r["abs_lambda_ci"])]
    if len(tied) == 1:
        return {
            "winner": best.get("experiment"),
            "reason": "unique_min_abs_lambda",
            "tied": [best.get("experiment")],
            "ranking_unstable": False,
        }
    with_ll = [r for r in tied if np.isfinite(r.get(ll_key, float("nan")))]
    if not with_ll:
        return {
            "winner": best.get("experiment"),
            "reason": "ci_overlap_no_ll",
            "tied": [r.get("experiment") for r in tied],
            "ranking_unstable": False,
        }
    with_ll.sort(key=lambda r: -float(r[ll_key]))
    return {
        "winner": with_ll[0].get("experiment"),
        "reason": "ci_overlap_mean_ll_tiebreak",
        "tied": [r.get("experiment") for r in tied],
        "ranking_unstable": False,
    }


def pick_winner_across_betas(
    by_beta: dict[str, list[dict[str, Any]]],
    *,
    preflight: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Require the same recipe to win at every β magnitude."""
    if preflight is not None:
        gate = assert_preflight_ok(preflight)
        if not gate["ok"]:
            return {
                "winner": None,
                "ranking_unstable": True,
                "reason": gate["reason"],
                "per_beta": {},
                "preflight": gate,
            }
    per_beta: dict[str, Any] = {}
    winners: list[Optional[str]] = []
    for label in sorted(by_beta.keys()):
        decision = pick_winner(by_beta[label], preflight=None)
        per_beta[label] = decision
        winners.append(decision.get("winner"))
    if not winners or any(w is None for w in winners):
        return {
            "winner": None,
            "ranking_unstable": True,
            "reason": "missing_winner",
            "per_beta": per_beta,
        }
    if len(set(winners)) == 1:
        return {
            "winner": winners[0],
            "ranking_unstable": False,
            "reason": "stable_across_betas",
            "per_beta": per_beta,
        }
    return {
        "winner": None,
        "ranking_unstable": True,
        "reason": "disagreement_across_betas",
        "per_beta": per_beta,
    }


def assert_preflight_ok(preflight: dict[str, Any]) -> dict[str, Any]:
    q = preflight.get("q_vs_pcs") or {}
    anc = preflight.get("anchor_or") or {}
    q_ok = bool(q.get("pass")) or bool(q.get("pcs_dropped"))
    # If PCs dropped, residual check on reduced set must pass
    if q.get("pcs_dropped") and not q.get("pass_after_drop", True):
        q_ok = False
    anc_ok = bool(anc.get("pass"))
    if not q_ok:
        return {"ok": False, "reason": "preflight_q_pcs_failed", "q_vs_pcs": q, "anchor_or": anc}
    if not anc_ok:
        return {"ok": False, "reason": "preflight_anchor_failed", "q_vs_pcs": q, "anchor_or": anc}
    return {"ok": True, "reason": "preflight_ok", "q_vs_pcs": q, "anchor_or": anc}


# ---------------------------------------------------------------------------
# Global ancestry I/O + regression helpers
# ---------------------------------------------------------------------------


def resolve_id_col(columns: Sequence[str]) -> str:
    for cand in ("ID", "IID", "sample", "SAMPLE", "research_id"):
        if cand in columns:
            return cand
    raise ValueError(f"no sample id column in {list(columns)}")


def load_global_anc(path: Path) -> "Any":
    """Load FLARE ``*.global.anc(.gz)`` or a TSV with ID + ancestry columns."""
    pd = _pd()
    with open_text(path) as fh:
        first = fh.readline()
    delim = "\t" if "\t" in first else None  # None → whitespace
    df = pd.read_csv(path, sep=delim, engine="python" if delim is None else "c")
    # normalize colnames
    rename = {c: c.lower() for c in df.columns}
    df = df.rename(columns=rename)
    id_col = None
    for cand in ("sample", "id", "iid", "research_id"):
        if cand in df.columns:
            id_col = cand
            break
    if id_col is None:
        raise SystemExit(f"{path}: need SAMPLE/ID column; got {list(df.columns)}")
    df = df.rename(columns={id_col: "ID"})
    df["ID"] = df["ID"].astype(str)
    missing = [a for a in ANC_NAMES if a not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing ancestry columns {missing}")
    for a in ANC_NAMES:
        df[a] = _pd().to_numeric(df[a], errors="coerce")
    return df[["ID", *ANC_NAMES]].drop_duplicates("ID")


def parse_covariates(raw: str) -> list[str]:
    return [c.strip() for c in raw.split(",") if c.strip()]


def pc_names(covariates: Sequence[str]) -> list[str]:
    return [c for c in covariates if re.fullmatch(r"PC\d+", c, flags=re.IGNORECASE) or c.upper().startswith("PC")]


def ols_r2(y: np.ndarray, X: np.ndarray) -> tuple[float, float]:
    """Return (R², residual SD) for y ~ 1 + X."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y = y[mask]
    X = X[mask]
    if y.size < 3:
        return float("nan"), float("nan")
    design = np.column_stack([np.ones(y.size), X])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(r2), float(np.sqrt(ss_res / max(y.size - design.shape[1], 1)))


def multivariate_r2(Q: np.ndarray, X: np.ndarray) -> dict[str, Any]:
    """Frobenius residual fraction + per-column R² for Q (n×k) ~ X."""
    Q = np.asarray(Q, dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    mask = np.all(np.isfinite(Q), axis=1) & np.all(np.isfinite(X), axis=1)
    Q = Q[mask]
    X = X[mask]
    per: dict[str, float] = {}
    resid_cols = []
    for j in range(Q.shape[1]):
        r2, _ = ols_r2(Q[:, j], X)
        per[str(j)] = r2
        # residual vector
        design = np.column_stack([np.ones(Q.shape[0]), X])
        beta, *_ = np.linalg.lstsq(design, Q[:, j], rcond=None)
        resid_cols.append(Q[:, j] - design @ beta)
    R = np.column_stack(resid_cols) if resid_cols else np.zeros_like(Q)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    ss_tot = float(np.sum(Qc**2))
    ss_res = float(np.sum(R**2))
    multi = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {
        "multivariate_r2": float(multi),
        "per_component_r2": per,
        "residual_frobenius": ss_res,
        "total_frobenius": ss_tot,
        "n": int(Q.shape[0]),
    }


def univariate_logit(y: np.ndarray, x: np.ndarray) -> dict[str, float]:
    """IRLS logistic regression y ~ 1 + x; returns slope, SE, Wald p."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(y) & np.isfinite(x)
    y = y[mask]
    x = x[mask]
    if y.size < 20 or y.min() == y.max():
        return {"beta": float("nan"), "se": float("nan"), "pval": float("nan"), "n": float(y.size)}
    X = np.column_stack([np.ones(y.size), x])
    beta = np.zeros(2)
    for _ in range(80):
        eta = np.clip(X @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(p * (1.0 - p), 1e-8, None)
        z = eta + (y - p) / W
        XtW = X.T * W
        try:
            beta_new = np.linalg.solve(XtW @ X, XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(beta_new - beta)) < 1e-9:
            beta = beta_new
            break
        beta = beta_new
    eta = np.clip(X @ beta, -30, 30)
    p = 1.0 / (1.0 + np.exp(-eta))
    W = np.clip(p * (1.0 - p), 1e-8, None)
    try:
        cov = np.linalg.inv((X.T * W) @ X)
        se = float(np.sqrt(max(cov[1, 1], 0.0)))
    except np.linalg.LinAlgError:
        se = float("nan")
    if not np.isfinite(se) or se <= 0:
        pval = float("nan")
    else:
        z = abs(float(beta[1]) / se)
        pval = float(math.erfc(z / math.sqrt(2.0)))
    return {"beta": float(beta[1]), "se": se, "pval": pval, "n": float(y.size)}


def merge_pheno_q(
    pheno_cov: Path,
    global_anc: Path,
    *,
    id_col: Optional[str] = None,
) -> tuple[Any, str]:
    pd = _pd()
    pheno = pd.read_csv(pheno_cov, sep="\t", dtype=str)
    idc = id_col or resolve_id_col(list(pheno.columns))
    pheno[idc] = pheno[idc].astype(str)
    q = load_global_anc(global_anc)
    merged = pheno.merge(q, left_on=idc, right_on="ID", how="inner", suffixes=("", "_q"))
    if merged.empty:
        raise SystemExit("no overlapping IDs between pheno_cov and global_anc")
    return merged, idc


def preflight_q_vs_pcs(
    pheno_cov: Path,
    global_anc: Path,
    covariates: Sequence[str],
    *,
    r2_max: float = DEFAULT_R2_MAX,
    primary_anc: str = "afr",
) -> dict[str, Any]:
    """Regress fixed q on the null PC set; drop PCs if R² too high."""
    pd = _pd()
    merged, idc = merge_pheno_q(pheno_cov, global_anc)
    pcs = [c for c in pc_names(list(covariates)) if c in merged.columns]
    non_pc = [c for c in covariates if c not in pcs]
    Q = merged[list(ANC_NAMES)].to_numpy(dtype=float)

    def _eval(pc_list: list[str]) -> dict[str, Any]:
        if not pc_list:
            return {
                "multivariate_r2": 0.0,
                "per_component_r2": {a: 0.0 for a in ANC_NAMES},
                "primary_r2": 0.0,
                "primary_resid_sd": float(np.nanstd(merged[primary_anc].astype(float))),
                "n_pcs": 0,
                "pc_columns": [],
            }
        X = merged[pc_list].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        multi = multivariate_r2(Q, X)
        # map per-component to anc names
        per_named = {
            ANC_NAMES[int(k)]: float(v) for k, v in multi["per_component_r2"].items() if int(k) < len(ANC_NAMES)
        }
        r2_p, resid_sd = ols_r2(
            merged[primary_anc].to_numpy(dtype=float),
            X,
        )
        return {
            "multivariate_r2": multi["multivariate_r2"],
            "per_component_r2": per_named,
            "primary_anc": primary_anc,
            "primary_r2": r2_p,
            "primary_resid_sd": resid_sd,
            "n_pcs": len(pc_list),
            "pc_columns": pc_list,
            "n": multi["n"],
        }

    with_pcs = _eval(pcs)
    passed = bool(with_pcs["multivariate_r2"] < r2_max)
    pcs_dropped = False
    after_drop = None
    null_covariates = list(covariates)
    if not passed:
        pcs_dropped = True
        null_covariates = list(non_pc)
        after_drop = _eval([])  # no PCs → R²=0 by construction
        # Remediating by drop always yields residual q
        passed_after = True
    else:
        passed_after = True

    report = {
        "pass": passed,
        "pcs_dropped": pcs_dropped,
        "pass_after_drop": passed_after if pcs_dropped else None,
        "r2_max": r2_max,
        "with_pcs": with_pcs,
        "after_drop": after_drop,
        "null_covariates_used": null_covariates,
        "requested_covariates": list(covariates),
        "n_samples_merged": int(len(merged)),
        "id_col": idc,
    }
    # Effective pass for downstream: original pass OR successful PC drop
    report["ok_for_pilot"] = bool(passed or (pcs_dropped and passed_after))
    return report


def preflight_anchor_or(
    pheno_cov: Path,
    global_anc: Path,
    candidates: Sequence[str],
    *,
    primary_anc: str = "afr",
    abs_beta_min: float = DEFAULT_OR_ABS_MIN,
    pmax: float = DEFAULT_OR_PMAX,
) -> dict[str, Any]:
    """Require non-trivial logit(anchor ~ primary q) before locking β_mid."""
    pd = _pd()
    merged, idc = merge_pheno_q(pheno_cov, global_anc)
    q = pd.to_numeric(merged[primary_anc], errors="coerce").to_numpy(dtype=float)
    trials = []
    chosen = None
    for pheno in candidates:
        if pheno not in merged.columns:
            trials.append({"phenotype": pheno, "status": "missing_column"})
            continue
        y = pd.to_numeric(merged[pheno], errors="coerce").to_numpy(dtype=float)
        fit = univariate_logit(y, q)
        ok = (
            np.isfinite(fit["beta"])
            and abs(fit["beta"]) >= abs_beta_min
            and np.isfinite(fit["pval"])
            and fit["pval"] < pmax
        )
        case_rate = float(np.nanmean(y == 1)) if np.isfinite(y).any() else float("nan")
        row = {
            "phenotype": pheno,
            "status": "pass" if ok else "fail",
            "log_or": fit["beta"],
            "se": fit["se"],
            "pval": fit["pval"],
            "n": fit["n"],
            "case_rate": case_rate,
            "abs_beta_min": abs_beta_min,
            "pmax": pmax,
            "primary_anc": primary_anc,
        }
        trials.append(row)
        if ok and chosen is None:
            chosen = row

    return {
        "pass": chosen is not None,
        "chosen": chosen,
        "trials": trials,
        "primary_anc": primary_anc,
        "abs_beta_min": abs_beta_min,
        "pmax": pmax,
        "id_col": idc,
        "beta_mid": None if chosen is None else float(chosen["log_or"]),
        "anchor_phenotype": None if chosen is None else chosen["phenotype"],
        "case_rate": None if chosen is None else float(chosen["case_rate"]),
    }


def run_preflight(
    *,
    pheno_cov: Path,
    global_anc: Path,
    covariates: Sequence[str],
    anchor_candidates: Sequence[str],
    primary_anc: str = "afr",
    r2_max: float = DEFAULT_R2_MAX,
    abs_beta_min: float = DEFAULT_OR_ABS_MIN,
    pmax: float = DEFAULT_OR_PMAX,
) -> dict[str, Any]:
    q_report = preflight_q_vs_pcs(
        pheno_cov, global_anc, covariates, r2_max=r2_max, primary_anc=primary_anc
    )
    # For gating pick_winner: treat pass=True if ok_for_pilot (incl. PC drop)
    q_report["pass"] = bool(q_report.get("ok_for_pilot"))
    anc_report = preflight_anchor_or(
        pheno_cov,
        global_anc,
        anchor_candidates,
        primary_anc=primary_anc,
        abs_beta_min=abs_beta_min,
        pmax=pmax,
    )
    report = {
        "q_vs_pcs": q_report,
        "anchor_or": anc_report,
        "null_covariates_used": q_report["null_covariates_used"],
        "beta_mid": anc_report.get("beta_mid"),
        "anchor_phenotype": anc_report.get("anchor_phenotype"),
        "case_rate": anc_report.get("case_rate"),
        "primary_anc": primary_anc,
    }
    report["ok"] = assert_preflight_ok(report)["ok"]
    return report


def beta_grid_from_mid(beta_mid: float) -> dict[str, float]:
    return {
        "lo": 0.5 * float(beta_mid),
        "mid": float(beta_mid),
        "hi": 2.0 * float(beta_mid),
    }


def simulate_structured_pheno(
    pheno_cov: Path,
    global_anc: Path,
    out_path: Path,
    *,
    beta: float,
    primary_anc: str = "afr",
    family_col: str = "pedigree_family_id",
    batch_col: str = "GC",
    case_rate: float,
    phenotype_col: str = SIM_PHENO_COL,
    sigma_f: float = 0.5,
    sigma_eps: float = 1.0,
    gc_effect_sd: float = 0.3,
    seed: int = 1,
) -> dict[str, Any]:
    """Simulate binary phenotype from liability with family shared RE.

    L = β·q_primary + γ_GC + u_family + ε; threshold to ``case_rate``.
    Unrelateds (missing/unique family ids) get singleton family draws.
    """
    pd = _pd()
    rng = np.random.default_rng(seed)
    merged, idc = merge_pheno_q(pheno_cov, global_anc)
    n = len(merged)
    q = pd.to_numeric(merged[primary_anc], errors="coerce").fillna(0.0).to_numpy(dtype=float)

    if batch_col not in merged.columns:
        raise SystemExit(f"batch column {batch_col!r} missing from pheno_cov")
    batches = merged[batch_col].astype(str).fillna("NA")
    uniq_b = sorted(batches.unique())
    gamma_map = {b: float(rng.normal(0.0, gc_effect_sd)) for b in uniq_b}
    gamma = np.array([gamma_map[b] for b in batches], dtype=float)

    if family_col in merged.columns:
        fam_raw = merged[family_col]
        fam_ids = []
        for i, v in enumerate(fam_raw):
            s = "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v).strip()
            if s in ("", "nan", "None", "NA", "<NA>"):
                fam_ids.append(f"__singleton_{merged.iloc[i][idc]}")
            else:
                fam_ids.append(s)
    else:
        fam_ids = [f"__singleton_{sid}" for sid in merged[idc].astype(str)]

    uniq_f = sorted(set(fam_ids))
    u_map = {f: float(rng.normal(0.0, sigma_f)) for f in uniq_f}
    u = np.array([u_map[f] for f in fam_ids], dtype=float)
    eps = rng.normal(0.0, sigma_eps, size=n)
    liability = float(beta) * q + gamma + u + eps

    # threshold so P(L > t) ≈ case_rate
    case_rate = float(np.clip(case_rate, 1e-4, 1.0 - 1e-4))
    t = float(np.quantile(liability, 1.0 - case_rate))
    y = (liability > t).astype(int)

    # write full pheno_cov clone with sim column (all original rows; missing sim=NA)
    pheno = pd.read_csv(pheno_cov, sep="\t", dtype=str)
    id_full = resolve_id_col(list(pheno.columns))
    pheno[id_full] = pheno[id_full].astype(str)
    sim_map = dict(zip(merged[idc].astype(str), y.astype(int)))
    pheno[phenotype_col] = pheno[id_full].map(lambda s: sim_map.get(s, ""))
    # drop rows without simulated phenotype for fit_null complete-cases
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pheno.to_csv(out_path, sep="\t", index=False)

    # Family structure diagnostic: shared u_f raises within-family covariance.
    multi = [f for f in uniq_f if fam_ids.count(f) >= 2]
    fam_means = []
    for f in multi:
        idx = [i for i, ff in enumerate(fam_ids) if ff == f]
        fam_means.append(float(np.mean(liability[idx])))
    if fam_means and np.var(liability) > 0:
        between_var = float(np.var(fam_means))
        icc_proxy = between_var / float(np.var(liability))
    else:
        icc_proxy = float("nan")
    c = liability - liability.mean()
    within_prod = []
    for f in multi:
        idx = [i for i, ff in enumerate(fam_ids) if ff == f]
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                within_prod.append(float(c[idx[a]] * c[idx[b]]))
    mean_within_prod = float(np.mean(within_prod)) if within_prod else float("nan")

    return {
        "path": str(out_path),
        "phenotype_col": phenotype_col,
        "beta": float(beta),
        "primary_anc": primary_anc,
        "case_rate_target": case_rate,
        "case_rate_obs": float(y.mean()),
        "n": int(n),
        "n_families_shared": int(sum(1 for f in uniq_f if not f.startswith("__singleton_"))),
        "sigma_f": sigma_f,
        "sigma_eps": sigma_eps,
        "family_icc_proxy": icc_proxy,
        "mean_within_family_centered_prod": mean_within_prod,
        "threshold": t,
    }


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, text=True, **kwargs)


def slice_anc_vcf(anc_vcf: str, region: str, out_vcf: Path) -> Path:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    run(["bcftools", "view", "-r", region, "-Oz", "-o", str(out_vcf), anc_vcf])
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
    ordered_dir = output_dir / "ordered"
    ordered_dir.mkdir(parents=True, exist_ok=True)
    dosage_files: list[Path] = []
    for p in sorted(output_dir.glob("*.dosage.txt.gz")):
        if ".collapse." in p.name or ".split." in p.name:
            continue
        m = re.search(r"(?<![.])\.anc(\d+)\.dosage\.txt\.gz$", p.name)
        if not m:
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
    for cand in (p_col, "p", "Pvalue", "pvalue", "P_CCT", "p_CCT"):
        if cand in df.columns:
            return _pd().to_numeric(df[cand], errors="coerce").to_numpy(dtype=float)
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
        "n_rep": len(list(lambdas)),
        "mean_lambda": float(np.nanmean(lambdas)) if len(lambdas) else float("nan"),
        "abs_lambda_dev": ci["mean"],
        "abs_lambda_ci": ci,
        "n_tested_sites": n_tested_sites,
        "min_tested_sites": min_tested_sites,
        "unstable_lambda": bool(unstable),
    }


def _add_shared_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--pheno-cov", type=Path, required=True)
    p.add_argument("--global-anc", type=Path, required=True)
    p.add_argument("--primary-anc", default="afr")
    p.add_argument("--family-col", default="pedigree_family_id")
    p.add_argument("--batch-col", default="GC")
    p.add_argument("--seed", type=int, default=1)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd")

    # ---- preflight ----
    pf = sub.add_parser("preflight", help="Run q↔PC and anchor-OR preflights")
    _add_shared_args(pf)
    pf.add_argument("--covariates", required=True, help="comma-separated null covariates")
    pf.add_argument(
        "--anchor-candidates",
        required=True,
        help="comma-separated binary phenotype columns to try as anchor",
    )
    pf.add_argument("--r2-max", type=float, default=DEFAULT_R2_MAX)
    pf.add_argument("--abs-beta-min", type=float, default=DEFAULT_OR_ABS_MIN)
    pf.add_argument("--pmax", type=float, default=DEFAULT_OR_PMAX)
    pf.add_argument("--out", type=Path, required=True)

    # ---- simulate-only (debug) ----
    sim = sub.add_parser("simulate", help="Write one simulated pheno_cov and exit")
    _add_shared_args(sim)
    sim.add_argument("--beta", type=float, required=True)
    sim.add_argument("--case-rate", type=float, required=True)
    sim.add_argument("--out", type=Path, required=True)
    sim.add_argument("--sigma-f", type=float, default=0.5)
    sim.add_argument("--sigma-eps", type=float, default=1.0)
    sim.add_argument("--gc-effect-sd", type=float, default=0.3)

    # ---- declare winner ----
    win = sub.add_parser("declare-winner", help="Pick winner across recipes / β grid")
    win.add_argument(
        "--scores",
        type=Path,
        required=True,
        help="JSON: either list of rows, or {beta_label: [rows]} for multi-β",
    )
    win.add_argument("--preflight", type=Path, default=None)
    win.add_argument("--out-dir", type=Path, required=True)

    # ---- aggregate-only ----
    agg = sub.add_parser("aggregate-only", help="Aggregate λ list without Tractor")
    agg.add_argument("--lambdas-json", type=Path, required=True)
    agg.add_argument("--out-dir", type=Path, required=True)
    agg.add_argument("--min-tested-sites", type=int, default=100)
    agg.add_argument("--seed", type=int, default=1)
    agg.add_argument("--experiment", default="")

    # ---- full run ----
    run_p = sub.add_parser("run", help="Extract tracts + simulated-pheno Tractor null-λ")
    _add_shared_args(run_p)
    run_p.add_argument("--anc-vcf", required=True)
    run_p.add_argument("--region", required=True)
    run_p.add_argument("--samples", required=True, type=Path)
    run_p.add_argument("--preflight", type=Path, required=True, help="null_lambda_preflight.json")
    run_p.add_argument(
        "--covariates",
        default="",
        help="optional override; default = preflight null_covariates_used",
    )
    run_p.add_argument("--grm-rds", required=True, type=Path)
    run_p.add_argument("--fit-null-r", type=Path, default=Path("scripts/fit_null.R"))
    run_p.add_argument("--extract-bin", default="extract-tracts-flare")
    run_p.add_argument("--score-bin", default="tractor-mix-score")
    run_p.add_argument("--out-dir", required=True, type=Path)
    run_p.add_argument("--n-rep", type=int, default=20)
    run_p.add_argument("--ac-threshold", type=int, default=20)
    run_p.add_argument("--min-tested-sites", type=int, default=100)
    run_p.add_argument("--num-ancs", type=int, default=5)
    run_p.add_argument("--threads", type=int, default=8)
    run_p.add_argument("--experiment", default="")
    run_p.add_argument("--sigma-f", type=float, default=0.5)
    run_p.add_argument("--sigma-eps", type=float, default=1.0)
    run_p.add_argument("--gc-effect-sd", type=float, default=0.3)
    run_p.add_argument(
        "--beta-scales",
        default="",
        help="optional comma lo,mid,hi overrides; default from preflight β_mid bracket",
    )
    run_p.add_argument("--dry-run", action="store_true")

    # Back-compat: legacy flat flags → error with migration hint
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] not in (
        "preflight",
        "simulate",
        "declare-winner",
        "aggregate-only",
        "run",
        "-h",
        "--help",
    ):
        if any(a.startswith("--") for a in argv):
            raise SystemExit(
                "flare_lai_null_lambda.py now uses subcommands: "
                "preflight | simulate | run | declare-winner | aggregate-only. "
                "Unstratified --phenotype shuffle has been removed."
            )

    args = parser.parse_args(argv)
    if not args.cmd:
        parser.print_help()
        return 2

    if args.cmd == "preflight":
        report = run_preflight(
            pheno_cov=args.pheno_cov,
            global_anc=args.global_anc,
            covariates=parse_covariates(args.covariates),
            anchor_candidates=parse_covariates(args.anchor_candidates),
            primary_anc=args.primary_anc,
            r2_max=args.r2_max,
            abs_beta_min=args.abs_beta_min,
            pmax=args.pmax,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"ok": report["ok"], "out": str(args.out)}, indent=2))
        if not report["ok"]:
            print("PREFLIGHT FAILED", file=sys.stderr)
            return 1
        return 0

    if args.cmd == "simulate":
        meta = simulate_structured_pheno(
            args.pheno_cov,
            args.global_anc,
            args.out,
            beta=args.beta,
            primary_anc=args.primary_anc,
            family_col=args.family_col,
            batch_col=args.batch_col,
            case_rate=args.case_rate,
            sigma_f=args.sigma_f,
            sigma_eps=args.sigma_eps,
            gc_effect_sd=args.gc_effect_sd,
            seed=args.seed,
        )
        print(json.dumps(meta, indent=2))
        return 0

    if args.cmd == "declare-winner":
        payload = json.loads(args.scores.read_text())
        preflight = None
        if args.preflight is not None:
            preflight = json.loads(args.preflight.read_text())
        args.out_dir.mkdir(parents=True, exist_ok=True)
        if isinstance(payload, dict) and all(isinstance(v, list) for v in payload.values()):
            decision = pick_winner_across_betas(payload, preflight=preflight)
        else:
            decision = pick_winner(payload, preflight=preflight)
        out = args.out_dir / "winner.json"
        out.write_text(json.dumps(decision, indent=2) + "\n")
        print(json.dumps(decision, indent=2))
        print("wrote", out)
        return 0 if not decision.get("ranking_unstable") else 1

    if args.cmd == "aggregate-only":
        payload = json.loads(args.lambdas_json.read_text())
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

    # ---- run ----
    preflight = json.loads(args.preflight.read_text())
    gate = assert_preflight_ok(preflight)
    if not gate["ok"]:
        raise SystemExit(f"preflight not OK: {gate['reason']}")

    null_covs = parse_covariates(args.covariates) if args.covariates else list(
        preflight.get("null_covariates_used") or []
    )
    if not null_covs:
        raise SystemExit("no null covariates (pass --covariates or preflight null_covariates_used)")

    beta_mid = preflight.get("beta_mid")
    case_rate = preflight.get("case_rate")
    if beta_mid is None or case_rate is None:
        raise SystemExit("preflight missing beta_mid / case_rate")
    if args.beta_scales.strip():
        parts = [float(x) for x in args.beta_scales.split(",")]
        if len(parts) != 3:
            raise SystemExit("--beta-scales needs lo,mid,hi")
        betas = {"lo": parts[0], "mid": parts[1], "hi": parts[2]}
    else:
        betas = beta_grid_from_mid(float(beta_mid))

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    plan = {
        "experiment": args.experiment,
        "anc_vcf": args.anc_vcf,
        "region": args.region,
        "n_rep": args.n_rep,
        "seed": args.seed,
        "betas": betas,
        "null_covariates_used": null_covs,
        "preflight_ok": True,
        "phenotype": SIM_PHENO_COL,
        "steps": [
            "bcftools view -r REGION anc.vcf",
            "extract-tracts-flare",
            "for each β in {lo,mid,hi}: for k in 1..K simulate liability pheno → fit_null → score → λGC",
            "aggregate per β; multi-β winner requires stable ranking",
        ],
    }
    (out_dir / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    if args.dry_run:
        print("wrote", out_dir / "plan.json")
        return 0

    for bin_name in (args.extract_bin, args.score_bin):
        if shutil.which(bin_name) is None:
            raise SystemExit(
                f"{bin_name} not found on PATH; use --dry-run or aggregate-only for offline scoring"
            )

    sliced = slice_anc_vcf(args.anc_vcf, args.region, out_dir / "window.anc.vcf.gz")
    dosage_files = extract_tracts(
        sliced,
        args.samples,
        out_dir / "extract",
        num_ancs=args.num_ancs,
        extract_bin=args.extract_bin,
    )

    per_beta_summaries: dict[str, Any] = {}
    for beta_idx, (label, beta) in enumerate(betas.items()):
        lambdas: list[float] = []
        n_tested_sites = 0
        for k in range(args.n_rep):
            seed_k = args.seed + 1000 * beta_idx + k
            rep_dir = out_dir / f"beta_{label}" / f"rep_{k:03d}"
            rep_dir.mkdir(parents=True, exist_ok=True)
            meta = simulate_structured_pheno(
                args.pheno_cov,
                args.global_anc,
                rep_dir / "pheno_cov.tsv",
                beta=beta,
                primary_anc=args.primary_anc,
                family_col=args.family_col,
                batch_col=args.batch_col,
                case_rate=float(case_rate),
                sigma_f=args.sigma_f,
                sigma_eps=args.sigma_eps,
                gc_effect_sd=args.gc_effect_sd,
                seed=seed_k,
            )
            (rep_dir / "sim_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
            null_rds = rep_dir / "null.rds"
            null_export = rep_dir / "null_export"
            fit_null(
                fit_null_r=args.fit_null_r,
                pheno_cov=rep_dir / "pheno_cov.tsv",
                phenotype=SIM_PHENO_COL,
                covariates=",".join(null_covs),
                grm_rds=args.grm_rds,
                out_null_rds=null_rds,
                out_null_export=null_export,
            )
            score_tsv = rep_dir / "scores.tsv"
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
            (rep_dir / "lambda.json").write_text(
                json.dumps(
                    {"beta": label, "rep": k, "lambda_gc": lam, "n_p": int(np.isfinite(pvals).sum())},
                    indent=2,
                )
                + "\n"
            )
        summary = aggregate_perm_lambdas(
            lambdas,
            n_tested_sites=n_tested_sites,
            min_tested_sites=args.min_tested_sites,
            seed=args.seed + beta_idx,
        )
        summary.update(
            {
                "experiment": args.experiment,
                "region": args.region,
                "ac_threshold": args.ac_threshold,
                "phenotype": SIM_PHENO_COL,
                "beta_label": label,
                "beta": beta,
                "null_covariates_used": null_covs,
                "anchor_phenotype": preflight.get("anchor_phenotype"),
                "case_rate": case_rate,
            }
        )
        if summary["unstable_lambda"]:
            summary["unstable_note"] = (
                f"n_tested_sites={n_tested_sites} < min_tested_sites={args.min_tested_sites}; "
                "widen window or lower --ac-threshold before comparing λ"
            )
        per_beta_summaries[label] = summary
        beta_out = out_dir / f"null_lambda_score.{label}.json"
        beta_out.write_text(json.dumps(summary, indent=2) + "\n")

    # mid-β summary as the default score path for enrich_association_scores
    mid = per_beta_summaries.get("mid") or next(iter(per_beta_summaries.values()))
    combined = {
        **mid,
        "per_beta": per_beta_summaries,
        "betas": betas,
        "preflight": {
            "q_vs_pcs_pass": (preflight.get("q_vs_pcs") or {}).get("pass"),
            "pcs_dropped": (preflight.get("q_vs_pcs") or {}).get("pcs_dropped"),
            "anchor_phenotype": preflight.get("anchor_phenotype"),
            "beta_mid": beta_mid,
        },
        "null_covariates_used": null_covs,
    }
    out = out_dir / "null_lambda_score.json"
    out.write_text(json.dumps(combined, indent=2) + "\n")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
