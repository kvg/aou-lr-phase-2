"""Tests for qq_spa_screen.py (post-hoc QC figure, not a WDL task)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "qqs", Path(__file__).with_name("qq_spa_screen.py"))
qqs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qqs)


def _frame(n_null=5000, n_artifact=8, n_real=2, seed=0):
    """Null p-values, plus artifacts (SPA bailed) and real hits (SPA converged)."""
    rng = np.random.default_rng(seed)
    p = list(rng.uniform(size=n_null))
    tstat, var, isspa = [1.0] * n_null, [1.0] * n_null, [False] * n_null
    for _ in range(n_artifact):          # |Z| far above cutoff, SPA not converged
        p.append(1e-40); tstat.append(30.0); var.append(1.0); isspa.append(False)
    for _ in range(n_real):              # |Z| above cutoff, SPA converged
        p.append(1e-9); tstat.append(6.0); var.append(1.0); isspa.append(True)
    df = pd.DataFrame({"P_cct_admixed": p, "Tstat_ancALL": tstat,
                       "var_ancALL": var, "Is.SPA_ancALL": isspa})
    for k in range(1, 6):
        df[f"Tstat_anc{k}"], df[f"var_anc{k}"], df[f"Is.SPA_anc{k}"] = 0.1, 1.0, True
    return df


def test_screen_removes_artifacts_and_keeps_real_hits(tmp_path):
    df = _frame()
    s = qqs.build(df, "P_cct_admixed", 5, 2.0, 5e-8, tmp_path / "q.png")
    assert s["n_spa_fallback"] == 8
    assert s["n_sig_as_reported"] == 10      # 8 artifacts + 2 real
    assert s["n_sig_screened"] == 2          # only the SPA-converged hits survive
    assert (tmp_path / "q.png").exists()


def test_bulk_null_rows_are_not_flagged():
    # |Z| = 1 is below the cutoff, so Is.SPA false is the normal state there
    df = _frame(n_artifact=0, n_real=0)
    mask, _, _ = qqs._fjta.spa_fallback_mask(df, 5, 2.0)
    assert int(mask.sum()) == 0


def test_lambda_is_essentially_unchanged_by_the_screen():
    # the screen removes a tail spike, not inflation/deflation of the bulk
    df = _frame()
    s = qqs.build(df, "P_cct_admixed", 5, 2.0, 5e-8, Path("/tmp/_qq_test.png"))
    assert abs(s["lambda_as_reported"] - s["lambda_screened"]) < 0.01


def test_tail_is_never_thinned():
    exp = np.linspace(0, 6, 50000)
    obs = np.linspace(0, 6, 50000)
    e, o = qqs.thin(exp, obs, keep_above=2.0, max_points=1000)
    assert o[o > 2.0].size == obs[obs > 2.0].size   # every tail point retained
    assert o.size < obs.size                        # bulk was thinned


def test_thin_retains_the_last_bulk_point():
    exp = np.arange(1000, dtype=float) / 1000
    obs = exp.copy()
    e, o = qqs.thin(exp, obs, keep_above=2.0, max_points=10)
    assert o[-1] == obs[-1]


def test_lambda_gc_of_uniform_is_about_one():
    rng = np.random.default_rng(1)
    assert abs(qqs.lambda_gc(rng.uniform(size=200000)) - 1.0) < 0.02
