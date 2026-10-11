#!/usr/bin/env python3
"""Four-arm (2x2) comparison of local ancestry x repeat-dosage coding.

Both arms come from FELIX step 2 on the same loci, the same Step 1 null model,
the same GRM and the same samples; the runs differ only in the DS{k} coding
written by scripts/write_admixed_dosage_vcf.py (``--coding ru`` vs
``--coding carrier``). Within each run, the pooled test is FELIX's own
ancestry-collapsed pass (``p.value_ancALL``, identical to ``P_hom_admixed``):
the score test on sum_k DS{k} with null model, GRM, variance ratio and marker
set held fixed. So the 2x2 is two runs, four columns:

    A pooled_biallelic   carrier coding, p.value_ancALL
    B la_biallelic       carrier coding, P_cct_admixed_c
    C pooled_dosage      ru coding,      p.value_ancALL
    D la_dosage          ru coding,      P_cct_admixed_c

Local-ancestry benefit is B vs A and D vs C; repeat-dosage benefit is C vs A
and D vs B. A is the conventional-GWAS baseline (global ancestry adjustment
through the Step 1 null only, no local ancestry, carrier coding).

Caveats recorded in the summary rather than hidden: the LA arms combine K
per-ancestry tests through CCT, so their p-values are not 1-df and the
reported lambda (computed by mapping p to a 1-df chi-square, the usual
convention for a p-value distribution) is a distributional summary, not a
genomic-control factor for a 1-df statistic. Significance-count differences
between arms are not a power comparison on their own; read them with the
per-locus table and the calibration panels.
"""

from __future__ import annotations

import argparse
import json
from math import erfc
from pathlib import Path

import numpy as np
import pandas as pd

ARMS = [
    ("pooled_biallelic", "biallelic", "pooled"),
    ("la_biallelic", "biallelic", "local-ancestry"),
    ("pooled_dosage", "repeat dosage", "pooled"),
    ("la_dosage", "repeat dosage", "local-ancestry"),
]

# First match wins. The pooled pass is unconditional, so p.value_c_ancALL and
# p.value_ancALL agree; P_hom_admixed is the same statistic reported again.
POOLED_P = ["p.value_ancALL", "P_hom_admixed_c", "P_hom_admixed"]
LA_P = ["P_cct_admixed_c", "P_cct_admixed"]
POOLED_BETA = ["BETA_ancALL", "BETA_c_ancALL"]
KEY = ["CHR", "POS", "MarkerID"]


def pick(df: pd.DataFrame, candidates: list[str], what: str, path: Path) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise SystemExit(f"{path}: no {what} column (looked for {', '.join(candidates)})")


def numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def lambda_from_p(p: np.ndarray) -> float:
    """Median 1-df chi-square from p / 0.4549364. No scipy in the FELIX image."""
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size < 10:
        return float("nan")

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


def anc_p_columns(df: pd.DataFrame, num_ancs: int) -> dict[str, str]:
    """Per-ancestry p-value columns, conditional variant preferred."""
    out = {}
    for k in range(1, num_ancs + 1):
        for c in (f"p.value_c_anc{k}", f"p.value_anc{k}"):
            if c in df.columns:
                out[f"anc{k}"] = c
                break
    return out


def load(path: Path, label: str, num_ancs: int) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in KEY if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing key column(s) {missing}")
    pooled_p = pick(df, POOLED_P, "pooled p-value", path)
    la_p = pick(df, LA_P, "local-ancestry p-value", path)
    beta = next((c for c in POOLED_BETA if c in df.columns), None)
    anc_cols = anc_p_columns(df, num_ancs)

    out = df[KEY].copy()
    out["MarkerID"] = out["MarkerID"].astype(str)
    out[f"P_pooled_{label}"] = numeric(df[pooled_p])
    out[f"P_la_{label}"] = numeric(df[la_p])
    if beta:
        out[f"BETA_pooled_{label}"] = numeric(df[beta])
    for name, col in anc_cols.items():
        out[f"P_{name}_{label}"] = numeric(df[col])
    if anc_cols:
        out[f"P_minanc_{label}"] = out[[f"P_{n}_{label}" for n in anc_cols]].min(axis=1)
    # Carrier haplotype count: an allele count only under the carrier coding.
    if label == "biallelic" and "AC_Allele2_ancALL" in df.columns:
        out["n_carrier_hap"] = numeric(df["AC_Allele2_ancALL"])
    prov = {
        "file": str(path), "n_rows": int(len(df)),
        "pooled_p_column": pooled_p, "la_p_column": la_p,
        "beta_column": beta, "per_ancestry_p_columns": anc_cols,
    }
    return out, prov


def transitions(sig_from: pd.Series, sig_to: pd.Series) -> dict[str, int]:
    return {
        "gained": int((~sig_from & sig_to).sum()),
        "lost": int((sig_from & ~sig_to).sum()),
        "both": int((sig_from & sig_to).sum()),
        "neither": int((~sig_from & ~sig_to).sum()),
    }


def markdown(summary: dict, p_threshold: float) -> str:
    lines = [
        f"### Four-arm comparison — {summary['phenotype']}",
        "",
        f"{summary['n_loci_shared']:,} loci tested in both runs "
        f"(threshold P < {p_threshold:g}).",
        "",
        "| Arm | Coding | Ancestry | Loci tested | lambda | Significant |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for arm, coding, ancestry in ARMS:
        s = summary["arms"][arm]
        lam = "n/a" if s["lambda"] != s["lambda"] else f"{s['lambda']:.3f}"
        lines.append(
            f"| {arm} | {coding} | {ancestry} | {s['n_tested']:,} | {lam} | {s['n_significant']:,} |"
        )
    lines += [
        "",
        "| Contrast | What it isolates | Gained | Lost |",
        "| --- | --- | ---: | ---: |",
    ]
    labels = {
        "pooled_biallelic->la_biallelic": "local ancestry, biallelic coding",
        "pooled_dosage->la_dosage": "local ancestry, repeat dosage",
        "pooled_biallelic->pooled_dosage": "repeat dosage, pooled",
        "la_biallelic->la_dosage": "repeat dosage, local ancestry",
        "pooled_biallelic->la_dosage": "both, vs the baseline",
    }
    for key, what in labels.items():
        t = summary["contrasts"][key]
        lines.append(f"| {key.replace('->', ' -> ')} | {what} | {t['gained']} | {t['lost']} |")
    lines += ["", f"Baseline arm: pooled_biallelic. {summary['notes']['la_p_not_1df']}", ""]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dosage-tsv", type=Path, required=True, help="FELIX results, --coding ru run")
    ap.add_argument("--biallelic-tsv", type=Path, required=True, help="FELIX results, --coding carrier run")
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--out-summary-json", type=Path, required=True)
    ap.add_argument("--out-md", type=Path, default=None)
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    ap.add_argument("--num-ancs", type=int, default=5)
    args = ap.parse_args()

    dos, dos_prov = load(args.dosage_tsv, "dosage", args.num_ancs)
    bia, bia_prov = load(args.biallelic_tsv, "biallelic", args.num_ancs)

    merged = bia.merge(dos, on=KEY, how="inner", validate="one_to_one")
    only_bia = len(bia) - len(merged)
    only_dos = len(dos) - len(merged)

    arm_p = {
        "pooled_biallelic": "P_pooled_biallelic",
        "la_biallelic": "P_la_biallelic",
        "pooled_dosage": "P_pooled_dosage",
        "la_dosage": "P_la_dosage",
    }
    sig = {arm: merged[col] < args.p_threshold for arm, col in arm_p.items()}

    summary = {
        "phenotype": args.phenotype,
        "p_threshold": args.p_threshold,
        "n_loci_shared": int(len(merged)),
        "n_loci_biallelic_only": int(only_bia),
        "n_loci_dosage_only": int(only_dos),
        "arms": {
            arm: {
                "coding": coding,
                "ancestry": ancestry,
                "p_column": arm_p[arm],
                "n_tested": int(merged[arm_p[arm]].notna().sum()),
                "lambda": lambda_from_p(merged[arm_p[arm]].to_numpy(dtype=float)),
                "n_significant": int(sig[arm].sum()),
            }
            for arm, coding, ancestry in ARMS
        },
        "contrasts": {
            f"{a}->{b}": transitions(sig[a], sig[b])
            for a, b in [
                ("pooled_biallelic", "la_biallelic"),
                ("pooled_dosage", "la_dosage"),
                ("pooled_biallelic", "pooled_dosage"),
                ("la_biallelic", "la_dosage"),
                ("pooled_biallelic", "la_dosage"),
            ]
        },
        "sources": {"dosage": dos_prov, "biallelic": bia_prov},
        "notes": {
            "la_p_not_1df": (
                "The local-ancestry arms combine K per-ancestry tests through CCT, so their "
                "p-values are not 1-df; lambda is a summary of the p-value distribution under a "
                "1-df mapping, not a genomic-control factor for a 1-df statistic."
            ),
            "marker_set": (
                "Both runs are written from the same locus VCF and marked on the same carrier "
                "haplotypes, but FELIX's carrier QC counts samples with a non-zero DS{k}: a "
                "sample whose gain and loss cancel on one ancestry is a carrier under the "
                "carrier coding and not under the repeat-unit coding, so a locus near --min-mac "
                "can pass in one arm only. n_loci_*_only counts those; the table and every count "
                "above are computed on the intersection."
            ),
        },
    }

    merged["min_p_any_arm"] = merged[list(arm_p.values())].min(axis=1)
    merged = merged.sort_values("min_p_any_arm")
    front = [c for c in KEY + ["n_carrier_hap"] + list(arm_p.values()) if c in merged.columns]
    merged = merged[front + [c for c in merged.columns if c not in front]]
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_tsv, sep="\t", index=False, na_rep="NA")
    args.out_summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    if args.out_md:
        args.out_md.write_text(markdown(summary, args.p_threshold))
    print(json.dumps({k: summary[k] for k in ("phenotype", "n_loci_shared", "arms", "contrasts")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
