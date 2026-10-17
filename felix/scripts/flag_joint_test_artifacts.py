#!/usr/bin/env python3
"""Flag FELIX joint-test artifacts in a merged result table, for manual review.

This is a POST-PROCESSING QC step. It is not part of FelixGenome.wdl, it is not
a WDL task input, and it changes nothing about how a run executes. It reads a
merged FELIX result TSV and writes the same rows back with QC flag columns
added. It never drops a row.

Why it exists
-------------
FELIX's joint ancestry statistic can return values no test on the data could
produce -- on a 54-case phenotype we saw P_het_admixed_c = 1.7e-18 where
P_het_admixed was 0.223, and 4.0e-140 on the repeat arm. The mechanism is
get_jointScore_pvalue (FELIX src/SKAT_funcs.cpp:1583): when the K x K
ancestry-stratified score covariance is numerically singular it falls back to a
Moore-Penrose pseudoinverse and returns a finite, extremely small p rather than
NA. The smallest ancestry panel is where this bites -- in our cohort anc2 holds
~5.3% of haplotypes, so a 54-case phenotype puts ~4 case haplotypes in it.

Flags (all are REVIEW flags; rows are kept either way)
------------------------------------------------------
effect_size_mismatch
    The hit's own pooled effect estimate does not support its p-value:
    2*Phi(-|BETA_ancALL/SE_ancALL|) > --p-threshold while the joint p is below
    it. p_implied reproduces p.value_ancALL to ~3 s.f., so the pooled test is
    internally coherent and the joint statistic is the one disagreeing.
    CAVEAT: a genuine ancestry-specific effect with opposite directions across
    ancestries cancels in the pooled test and will also raise this flag. That
    is precisely why this tool flags rather than excludes.

degenerate_covariance
    An ancestry was dropped from the fit (var_anc{k} missing), or the diagonal
    of the ancestry covariance is near-singular
    (min(var_anc*)/max(var_anc*) < --var-ratio). Targets the mechanism itself,
    so it does not penalise genuine ancestry-specific signal.
    CAVEAT: the default --var-ratio is NOT calibrated against a null
    distribution. Treat it as a starting point and tune it on a high-case
    phenotype before relying on it.

impossible_pvalue
    A reported p-value is exactly 0. Not representable as a test result; seen
    on the unconditional branch too.

negative_dosage_af
    AF_case_ancALL < 0. Informational for the repeat arms only: repeat dosage
    is REF-relative, so FELIX's sum(DS)/2N is not an allele frequency. Not
    evidence of an artifact on its own.

spa_fallback
    SPA was attempted on this test and silently fell back to the uncorrected
    normal approximation, so the reported p-value is not case-control-imbalance
    corrected. Fires when |Tstat|/sqrt(var) > --spa-cutoff (the threshold that
    makes FELIX attempt SPA at all) while Is.SPA is false, in any ancestry
    stratum or the pooled test.

    Mechanism: SAIGE_test.cpp:776 enters the SPA branch on
    StdStat > m_SPA_Cutoff; if the saddlepoint root-find does not converge, or
    the SPA p-value underflows to exactly 0 in non-log mode
    (SAIGE_test.cpp:908-915), t_isSPAConverge is set false and
    SAIGE_test.cpp:938-941 reports t_pval_noSPA instead -- the raw score-test
    p-value -- with no indication beyond the Is.SPA column. On a signed,
    REF-relative repeat dosage the score variance can collapse onto one or two
    outlier haplotypes, and the uncorrected normal p is then arbitrarily small:
    on a 63-case phenotype this produced P = 3.9e-135 (chi2 = 612) at a locus
    whose biallelic test gave P = 0.96.

    Is.SPA false on its own is NOT a problem -- it is the normal state for a
    test whose score never reached the cutoff, which is most of them. The
    conjunction with a large StdStat is what identifies an attempted-and-failed
    correction. Validated on the top_hits.tsv of three real phenotypes
    (BI_160 / BI_160_11 / BI_174, 50 rows per arm): the flag fires on 0 of the
    300 main-arm and ru_biallelic-arm rows (3 phenotypes x 2 arms x 50), and on
    48 of the 50 genome-wide hits in the 63-case repeat-dosage arm under
    P_cct_admixed (47 of 49 under P_cct_admixed_c). Under either column the
    survivors are the same two markers at chr9:8.98Mb, which the independent
    SNV arm also reports at P = 1.4e-09.

Usage
-----
    ./flag_joint_test_artifacts.py --results merged.felix.tsv.gz \
        --out-tsv merged.flagged.tsv.gz --out-summary-json qc.json

By default only rows below --p-threshold are examined (the review set); pass
--all-rows to flag every row.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P_JOINT_PREFERRED = ["P_cct_admixed", "P_cct_admixed_c"]
FLAG_COLS = ["effect_size_mismatch", "degenerate_covariance",
             "impossible_pvalue", "negative_dosage_af", "spa_fallback"]

# FELIX/SAIGE default SPAcutoff; run_felix_step2.R does not override it.
SPA_CUTOFF_DEFAULT = 2.0


def as_bool(s: pd.Series) -> pd.Series:
    """FELIX writes Is.SPA with std::boolalpha, so it arrives as either real
    bools (pandas sniffed the column) or the strings 'true'/'false'. NA stays
    NA: a missing Is.SPA is unknown, not a failure."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False})


def spa_fallback_mask(df: pd.DataFrame, num_ancs: int, spa_cutoff: float):
    """Tests where SPA was attempted but the uncorrected p-value was reported.

    Returns (mask, n_strata, max_std_stat) over the per-ancestry tests and the
    pooled ancALL test.
    """
    strata = [f"anc{k}" for k in range(1, num_ancs + 1)] + ["ancALL"]
    mask = pd.Series(False, index=df.index)
    n_strata = pd.Series(0, index=df.index, dtype=int)
    max_std = pd.Series(np.nan, index=df.index, dtype=float)
    for a in strata:
        tc, vc, sc = f"Tstat_{a}", f"var_{a}", f"Is.SPA_{a}"
        if not {tc, vc, sc} <= set(df.columns):
            continue
        t = pd.to_numeric(df[tc], errors="coerce")
        v = pd.to_numeric(df[vc], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            std = (t.abs() / np.sqrt(v.where(v > 0)))
        hit = ((as_bool(df[sc]) == False) & (std > spa_cutoff)).fillna(False)  # noqa: E712
        mask |= hit
        n_strata += hit.astype(int)
        max_std = pd.concat([max_std, std.where(hit)], axis=1).max(axis=1)
    return mask, n_strata, max_std


def pick_joint_col(df: pd.DataFrame) -> str:
    for c in P_JOINT_PREFERRED:
        if c in df.columns:
            return c
    raise SystemExit(f"no joint p-value column found; looked for {P_JOINT_PREFERRED}")


def norm_sf(z: np.ndarray) -> np.ndarray:
    """Two-sided normal tail probability without a scipy dependency."""
    return np.array([math.erfc(abs(v) / math.sqrt(2.0)) if np.isfinite(v) else np.nan
                     for v in np.asarray(z, dtype=float)])


def flag(df: pd.DataFrame, pcol: str, p_threshold: float, var_ratio: float,
         num_ancs: int, all_rows: bool,
         spa_cutoff: float = SPA_CUTOFF_DEFAULT) -> pd.DataFrame:
    out = df.copy()
    p = pd.to_numeric(out[pcol], errors="coerce")
    considered = pd.Series(True, index=out.index) if all_rows else (p < p_threshold)

    # effect_size_mismatch
    if {"BETA_ancALL", "SE_ancALL"} <= set(out.columns):
        beta = pd.to_numeric(out["BETA_ancALL"], errors="coerce")
        se = pd.to_numeric(out["SE_ancALL"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            out["p_implied_from_beta_se"] = norm_sf(beta / se)
        out["effect_size_mismatch"] = considered & (out["p_implied_from_beta_se"] > p_threshold)
    else:
        out["p_implied_from_beta_se"] = np.nan
        out["effect_size_mismatch"] = False

    # degenerate_covariance
    vcols = [f"var_anc{k}" for k in range(1, num_ancs + 1) if f"var_anc{k}" in out.columns]
    if vcols:
        V = out[vcols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            vmin, vmax = np.nanmin(V, axis=1), np.nanmax(V, axis=1)
            ratio = np.where(vmax > 0, vmin / vmax, np.nan)
        out["n_anc_dropped"] = np.isnan(V).sum(axis=1)
        out["var_ratio_anc"] = ratio
        out["degenerate_covariance"] = considered & (
            (out["n_anc_dropped"] > 0) | (pd.Series(ratio, index=out.index) < var_ratio))
    else:
        out["n_anc_dropped"] = np.nan
        out["var_ratio_anc"] = np.nan
        out["degenerate_covariance"] = False

    # impossible_pvalue — checked on every p-value column present
    pcols = [c for c in out.columns
             if c.startswith(("P_cct", "P_het", "P_hom", "p.value")) and "NA" not in c]
    zero = pd.Series(False, index=out.index)
    for c in pcols:
        zero |= (pd.to_numeric(out[c], errors="coerce") == 0)
    out["impossible_pvalue"] = zero          # never gated: an exact 0 is wrong anywhere

    # negative_dosage_af
    out["negative_dosage_af"] = (
        pd.to_numeric(out["AF_case_ancALL"], errors="coerce") < 0
        if "AF_case_ancALL" in out.columns else False)

    # spa_fallback
    spa_mask, spa_n, spa_std = spa_fallback_mask(out, num_ancs, spa_cutoff)
    out["n_strata_spa_fallback"] = spa_n
    out["max_std_stat_spa_fallback"] = spa_std
    out["spa_fallback"] = considered & spa_mask

    for c in FLAG_COLS:
        out[c] = out[c].fillna(False).astype(bool)
    out["qc_flags"] = [",".join(f for f in FLAG_COLS if r[f]) or "." for _, r in out[FLAG_COLS].iterrows()]
    out["qc_n_flags"] = out[FLAG_COLS].sum(axis=1).astype(int)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, required=True, help="merged FELIX result TSV (.tsv or .tsv.gz)")
    ap.add_argument("--out-tsv", type=Path, required=True, help="same rows, QC columns appended")
    ap.add_argument("--out-summary-json", type=Path, default=None)
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    ap.add_argument("--var-ratio", type=float, default=1e-4,
                    help="UNCALIBRATED default; tune on a high-case phenotype")
    ap.add_argument("--num-ancs", type=int, default=5)
    ap.add_argument("--spa-cutoff", type=float, default=SPA_CUTOFF_DEFAULT,
                    help="FELIX's SPAcutoff; a test above this had SPA attempted")
    ap.add_argument("--all-rows", action="store_true",
                    help="flag every row, not just those below --p-threshold")
    a = ap.parse_args()

    df = pd.read_csv(a.results, sep="\t", low_memory=False)
    pcol = pick_joint_col(df)
    out = flag(df, pcol, a.p_threshold, a.var_ratio, a.num_ancs, a.all_rows,
               spa_cutoff=a.spa_cutoff)

    n_hits = int((pd.to_numeric(out[pcol], errors="coerce") < a.p_threshold).sum())
    flagged = {c: int(out[c].sum()) for c in FLAG_COLS}
    summary = dict(results=str(a.results), joint_p_col=pcol, p_threshold=a.p_threshold,
                   var_ratio=a.var_ratio, spa_cutoff=a.spa_cutoff,
                   n_rows=int(len(out)), n_hits=n_hits,
                   flagged=flagged,
                   n_hits_any_flag=int((out["qc_n_flags"].gt(0) &
                                        (pd.to_numeric(out[pcol], errors="coerce") < a.p_threshold)).sum()),
                   n_hits_unflagged=int((out["qc_n_flags"].eq(0) &
                                         (pd.to_numeric(out[pcol], errors="coerce") < a.p_threshold)).sum()))

    out.to_csv(a.out_tsv, sep="\t", index=False)
    if a.out_summary_json:
        a.out_summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {a.out_tsv} ({len(out)} rows, no rows dropped)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
