#!/usr/bin/env python3
"""Simulation harness for repeat-unit dosage vs collapse/split encodings.

Generates synthetic per-sample repeat copy numbers under three architectures:
  length_additive — effect scales with total copy number C
  single_allele   — effect only when ALT allele is present (collapse-like)
  threshold       — effect when C exceeds a cutoff

For each replicate, fits a simple score test on dosage / collapse / split
encodings and reports empirical type I error and power.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ARCHITECTURES = ("length_additive", "single_allele", "threshold")
ENCODINGS = ("dosage", "collapse", "split")


@dataclass
class SimConfig:
    n_samples: int
    n_reps: int
    effect: float
    null_effect: float
    period: int
    cn_ref: float
    alt_ru: int
    threshold: int
    seed: int


def _erfcinv_scalar(p: float) -> float:
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if math.erfc(mid) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def p_from_z(z: float) -> float:
    if not np.isfinite(z):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2.0))


def encode_dosage(cn: np.ndarray, cn_ref: float, alt_ru: int) -> np.ndarray:
    """Copy-number dosage: REF = cn_ref, ALT haplotype adds signed RU."""
    return cn.astype(float)


def encode_collapse(cn: np.ndarray, cn_ref: float, alt_ru: int) -> np.ndarray:
    alt_c = cn_ref + alt_ru
    return (cn > cn_ref).astype(float)


def encode_split(cn: np.ndarray, cn_ref: float, alt_ru: int) -> np.ndarray:
    alt_c = cn_ref + alt_ru
    return (cn == alt_c).astype(float)


def phenotype_from_arch(
    g: np.ndarray,
    architecture: str,
    effect: float,
    cn_ref: float,
    alt_ru: int,
    threshold: int,
    rng: np.random.Generator,
) -> np.ndarray:
    if architecture == "length_additive":
        signal = effect * (g - cn_ref)
    elif architecture == "single_allele":
        signal = effect * (g > cn_ref).astype(float)
    elif architecture == "threshold":
        signal = effect * (g >= threshold).astype(float)
    else:
        raise ValueError(architecture)
    noise = rng.normal(0.0, 1.0, size=g.shape[0])
    return signal + noise


def score_test(y: np.ndarray, x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 20 or np.std(x) < 1e-8:
        return float("nan")
    x_c = x - x.mean()
    y_c = y - y.mean()
    beta = float(x_c @ y_c / (x_c @ x_c))
    resid = y_c - beta * x_c
    se = math.sqrt(float((resid @ resid) / max(1, x.size - 2)) / float(x_c @ x_c))
    if se <= 0:
        return float("nan")
    z = beta / se
    return p_from_z(z)


def simulate_copy_numbers(n: int, cn_ref: float, alt_ru: int, rng: np.random.Generator) -> np.ndarray:
    """Draw diploid copy numbers at one locus (hom/ref, het, hom/alt)."""
    geno = rng.integers(0, 3, size=n)
    alt_c = cn_ref + alt_ru
    out = np.where(geno == 0, cn_ref, np.where(geno == 1, 0.5 * (cn_ref + alt_c), alt_c))
    return out


def run_replicates(cfg: SimConfig, architecture: str, effect: float) -> pd.DataFrame:
    rng = np.random.default_rng(cfg.seed)
    rows = []
    for rep in range(cfg.n_reps):
        cn = simulate_copy_numbers(cfg.n_samples, cfg.cn_ref, cfg.alt_ru, rng)
        y = phenotype_from_arch(
            cn, architecture, effect, cfg.cn_ref, cfg.alt_ru, cfg.threshold, rng
        )
        encoders = {
            "dosage": encode_dosage,
            "collapse": encode_collapse,
            "split": encode_split,
        }
        for enc_name, enc_fn in encoders.items():
            x = enc_fn(cn, cfg.cn_ref, cfg.alt_ru)
            p = score_test(y, x)
            rows.append(
                {
                    "rep": rep,
                    "architecture": architecture,
                    "encoding": enc_name,
                    "effect": effect,
                    "p_value": p,
                }
            )
    return pd.DataFrame(rows)


def summarize_frame(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    rows = []
    for (arch, enc, effect), sub in df.groupby(["architecture", "encoding", "effect"]):
        p = sub["p_value"].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        rows.append(
            {
                "architecture": arch,
                "encoding": enc,
                "effect": effect,
                "n_reps": len(sub),
                "n_finite": p.size,
                "reject_rate": float((p < alpha).mean()) if p.size else float("nan"),
                "median_p": float(np.median(p)) if p.size else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-samples", type=int, default=2000)
    ap.add_argument("--n-reps", type=int, default=500)
    ap.add_argument("--effect", type=float, default=0.15, help="Alternative effect size")
    ap.add_argument("--null-effect", type=float, default=0.0)
    ap.add_argument("--cn-ref", type=float, default=10.0)
    ap.add_argument("--alt-ru", type=int, default=3)
    ap.add_argument("--threshold", type=int, default=15, help="For threshold architecture")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    cfg = SimConfig(
        n_samples=args.n_samples,
        n_reps=args.n_reps,
        effect=args.effect,
        null_effect=args.null_effect,
        period=3,
        cn_ref=args.cn_ref,
        alt_ru=args.alt_ru,
        threshold=args.threshold,
        seed=args.seed,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for arch in ARCHITECTURES:
        parts.append(run_replicates(cfg, arch, cfg.null_effect))
        parts.append(run_replicates(cfg, arch, cfg.effect))
    raw = pd.concat(parts, ignore_index=True)
    raw.to_csv(args.out_dir / "simulation_pvalues.tsv", sep="\t", index=False)

    summary = summarize_frame(raw, args.alpha)
    summary.to_csv(args.out_dir / "simulation_summary.tsv", sep="\t", index=False)

    # Type I error rows (null effect)
    t1 = summary.loc[summary["effect"] == cfg.null_effect].copy()
    t1.rename(columns={"reject_rate": "type1_error"}, inplace=True)
    t1.to_csv(args.out_dir / "type1_error.tsv", sep="\t", index=False)

    # Power rows (alternative effect)
    pw = summary.loc[summary["effect"] == cfg.effect].copy()
    pw.rename(columns={"reject_rate": "power"}, inplace=True)
    pw.to_csv(args.out_dir / "power.tsv", sep="\t", index=False)

    print(f"Wrote {args.out_dir / 'simulation_summary.tsv'} ({len(summary)} rows)")


if __name__ == "__main__":
    main()
