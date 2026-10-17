"""Tests for flag_joint_test_artifacts.py (post-hoc QC, not a WDL task)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "fjta", Path(__file__).with_name("flag_joint_test_artifacts.py"))
fjta = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fjta)

THR = 5e-8


def _row(**kw):
    base = dict(MarkerID="m", P_cct_admixed=0.5, BETA_ancALL=0.1, SE_ancALL=0.1,
                AF_case_ancALL=0.2, p_value_ancALL=0.5)
    base.update(kw)
    for k in range(1, 6):
        base.setdefault(f"var_anc{k}", 1.0)
    return base


def _flag(rows):
    df = pd.DataFrame(rows)
    return fjta.flag(df, "P_cct_admixed", THR, 1e-4, 5, all_rows=False)


def test_clean_row_unflagged():
    out = _flag([_row()])
    assert out["qc_n_flags"].iloc[0] == 0
    assert out["qc_flags"].iloc[0] == "."


def test_effect_size_mismatch_fires_only_on_hits():
    # tiny joint p, but BETA/SE implies p ~ 0.32 -> mismatch
    hit = _row(P_cct_admixed=1e-20, BETA_ancALL=0.1, SE_ancALL=0.1)
    # same inconsistency but not a hit -> not flagged (review set is hits only)
    nonhit = _row(P_cct_admixed=0.5, BETA_ancALL=0.1, SE_ancALL=0.1)
    out = _flag([hit, nonhit])
    assert bool(out["effect_size_mismatch"].iloc[0])
    assert not bool(out["effect_size_mismatch"].iloc[1])


def test_effect_size_consistent_hit_not_flagged():
    # z = 20 -> p ~ 5.5e-89, consistent with a genuine tiny joint p
    out = _flag([_row(P_cct_admixed=1e-20, BETA_ancALL=2.0, SE_ancALL=0.1)])
    assert not bool(out["effect_size_mismatch"].iloc[0])


def test_degenerate_covariance_dropped_ancestry():
    out = _flag([_row(P_cct_admixed=1e-20, BETA_ancALL=2.0, SE_ancALL=0.1, var_anc2=np.nan)])
    assert bool(out["degenerate_covariance"].iloc[0])
    assert out["n_anc_dropped"].iloc[0] == 1


def test_degenerate_covariance_var_ratio():
    out = _flag([_row(P_cct_admixed=1e-20, BETA_ancALL=2.0, SE_ancALL=0.1, var_anc2=1e-9)])
    assert bool(out["degenerate_covariance"].iloc[0])


def test_impossible_pvalue_not_gated_on_hits():
    # exact zero is wrong wherever it appears, including on a non-hit row
    out = _flag([_row(P_cct_admixed=0.0)])
    assert bool(out["impossible_pvalue"].iloc[0])


def test_negative_dosage_af_is_informational():
    out = _flag([_row(AF_case_ancALL=-0.01)])
    assert bool(out["negative_dosage_af"].iloc[0])
    assert out["qc_n_flags"].iloc[0] == 1


def _spa_row(std_stat, is_spa, p=1e-20, anc="ancALL"):
    """A row whose Tstat/var give the requested |Tstat|/sqrt(var)."""
    r = _row(P_cct_admixed=p)
    r[f"var_{anc}"] = 1.0
    r[f"Tstat_{anc}"] = float(std_stat)
    r[f"Is.SPA_{anc}"] = is_spa
    return r


def test_spa_fallback_fires_when_spa_attempted_and_failed():
    # |Z| = 24.7 is far above the cutoff, so SPA was attempted; false means it
    # bailed and the uncorrected normal p was reported.
    out = _flag([_spa_row(24.7, False)])
    assert bool(out["spa_fallback"].iloc[0])
    assert out["n_strata_spa_fallback"].iloc[0] == 1
    assert out["max_std_stat_spa_fallback"].iloc[0] == 24.7


def test_spa_not_attempted_is_not_a_fallback():
    # the normal state for most tests: score below cutoff, so Is.SPA is false
    # but no correction was ever needed
    out = _flag([_spa_row(0.8, False)])
    assert not bool(out["spa_fallback"].iloc[0])
    assert out["n_strata_spa_fallback"].iloc[0] == 0


def test_spa_converged_is_not_a_fallback():
    out = _flag([_spa_row(24.7, True)])
    assert not bool(out["spa_fallback"].iloc[0])


def test_spa_fallback_detected_in_an_ancestry_stratum():
    # headline p is a CCT over per-ancestry tests, so a single bad stratum
    # must raise the flag even when the pooled test is clean
    row = _spa_row(0.5, True)
    row.update({"var_anc3": 1.0, "Tstat_anc3": 30.0, "Is.SPA_anc3": False})
    out = _flag([row])
    assert bool(out["spa_fallback"].iloc[0])
    assert out["n_strata_spa_fallback"].iloc[0] == 1


def test_spa_fallback_accepts_boolalpha_strings():
    # FELIX writes std::boolalpha, so the column can arrive as text
    out = _flag([_spa_row(24.7, "false"), _spa_row(24.7, "true")])
    assert bool(out["spa_fallback"].iloc[0])
    assert not bool(out["spa_fallback"].iloc[1])


def test_spa_fallback_gated_on_the_review_set():
    out = _flag([_spa_row(24.7, False, p=0.5)])
    assert not bool(out["spa_fallback"].iloc[0])


def test_no_rows_are_dropped():
    rows = [_row(P_cct_admixed=p) for p in (1e-20, 0.5, 0.0, 1e-3)]
    out = _flag(rows)
    assert len(out) == len(rows)


def test_missing_optional_columns_do_not_crash():
    df = pd.DataFrame([{"MarkerID": "m", "P_cct_admixed": 1e-20}])
    out = fjta.flag(df, "P_cct_admixed", THR, 1e-4, 5, all_rows=False)
    assert len(out) == 1 and out["qc_n_flags"].iloc[0] == 0


def test_pick_joint_col_prefers_unconditional():
    df = pd.DataFrame([{"P_cct_admixed_c": 1e-9, "P_cct_admixed": 0.2}])
    assert fjta.pick_joint_col(df) == "P_cct_admixed"
