"""Tests for sweep_best_p.py."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location("sbp", Path(__file__).with_name("sweep_best_p.py"))
sbp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sbp)


def _top(rows):
    """rows: (marker, p, tstat, var, is_spa, carriers)."""
    return pd.DataFrame([{"MarkerID": m, "CHR": "chr9", "POS": 100 + i, "P_cct_admixed_c": p, "Tstat_ancALL": t, "var_ancALL": v,
                          "Is.SPA_ancALL": s, "N_case_het_ancALL": c, "N_case_hom_ancALL": 0} for i, (m, p, t, v, s, c) in enumerate(rows)])


def _bundle():
    main = _top([("fallback", 1e-15, 30.0, 1.0, False, 20),  # SPA attempted and failed
                 ("nocarrier", 1e-13, 5.0, 1.0, True, 0),
                 ("zero", 0.0, 5.0, 1.0, True, 9),
                 ("good", 1e-9, 6.0, 1.0, True, 20),
                 ("weak", 1e-5, 1.0, 1.0, False, 30)])  # Is.SPA false is normal below the cutoff
    ru = _top([("a", 1e-30, 30.0, 1.0, False, 12), ("b", 1e-20, 25.0, 1.0, False, 12), ("c", 1e-3, 1.0, 1.0, True, 0)])
    return {"arms": {"main": {"top_hits": main, "phenotype": "BI_174"}, "ru": {"top_hits": ru, "phenotype": "BI_174"}, "ru_biallelic": {}}}


def test_best_screened_row_skips_artifacts():
    r = sbp.best_p_report(_bundle()).set_index("arm")
    assert list(r.index) == ["main", "ru"]  # ru_biallelic has no top_hits
    m = r.loc["main"]
    assert m.MarkerID == "good" and m.best_p == 1e-9 and m.case_carriers == 20 and not m.censored
    assert m.n_rows == 5 and m.n_screened_in == 2  # good and weak survive


def test_all_rows_failing_is_censored_not_zero():
    ru = sbp.best_p_report(_bundle()).set_index("arm").loc["ru"]
    assert ru.censored and np.isnan(ru.best_p) and ru.n_screened_in == 0
    assert ru.p_bound == 1e-3  # the weakest P in the file


def test_sweep_reads_the_harvest_cache(tmp_path):
    for run in ("felix.prod.BI_174", "felix.prod.BI_160"):
        d = tmp_path / run
        (d / "main" / "qc" / "BI_174").mkdir(parents=True)
        (d / ".done").write_text("ok")
        _bundle()["arms"]["main"]["top_hits"].to_csv(d / "main" / "qc" / "BI_174" / "top_hits.tsv", sep="\t", index=False)
    (tmp_path / "felix.prod.partial").mkdir()  # no .done: ignored
    df = sbp.sweep_best_p(str(tmp_path))
    assert list(df.columns) == sbp.COLUMNS
    assert sorted(df.run.unique()) == ["felix.prod.BI_160", "felix.prod.BI_174"] and set(df.arm) == {"main"}
    assert (df.best_p == 1e-9).all()
