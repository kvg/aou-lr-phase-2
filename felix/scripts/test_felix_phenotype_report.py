"""Tests for felix_phenotype_report.py (post-hoc QC, not a WDL task).

Includes a drift guard: the notebook embeds this module via %%writefile, so
the two copies must stay byte-identical.
"""
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "fpr", Path(__file__).with_name("felix_phenotype_report.py"))
fpr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fpr)

NOTEBOOK = (Path(__file__).resolve().parents[2]
            / "notebooks" / "terra" / "felix_07_phenotype_review.ipynb")


# ---------------------------------------------------------------- drift guard

def test_notebook_embeds_this_exact_module():
    """The %%writefile cell must match felix/scripts/felix_phenotype_report.py.

    If this fails, the notebook and the committed module have diverged — the
    notebook is regenerated from the module, never edited in place.
    """
    assert NOTEBOOK.exists(), f"notebook not found at {NOTEBOOK}"
    nb = json.loads(NOTEBOOK.read_text())
    cells = [c for c in nb["cells"]
             if c["cell_type"] == "code"
             and "".join(c["source"]).startswith("%%writefile felix_phenotype_report.py")]
    assert len(cells) == 1, f"expected exactly one %%writefile cell, found {len(cells)}"
    embedded = "".join(cells[0]["source"]).split("\n", 1)[1]
    on_disk = Path(__file__).with_name("felix_phenotype_report.py").read_text()
    assert embedded == on_disk, "notebook copy has drifted from the module"


# ---------------------------------------------------------------- unit tests

def test_as_bool_handles_boolalpha_and_na():
    s = fpr.as_bool(pd.Series(["true", "FALSE", "True", None]))
    assert list(s[:3]) == [True, False, True]
    assert pd.isna(s.iloc[3])


def test_lambda_gc_of_uniform_is_about_one():
    rng = np.random.default_rng(0)
    assert abs(fpr.lambda_gc(rng.uniform(size=200_000)) - 1.0) < 0.02


def _frame(**over):
    row = {"MarkerID": "m", "N_case_het_ancALL": 10, "N_case_hom_ancALL": 2}
    for a in [f"anc{k}" for k in range(1, 6)] + ["ancALL"]:
        row |= {f"Tstat_{a}": 0.1, f"var_{a}": 1.0, f"Is.SPA_{a}": True,
                f"Tstat_c_{a}": 0.1, f"var_c_{a}": 1.0,
                f"p.value_c_{a}": 0.5, f"p.value.NA_c_{a}": 0.4}
    row.update(over)
    return pd.DataFrame([row])


def test_unconditioned_fallback_detected():
    m, which = fpr.spa_fallback(_frame(**{"Tstat_ancALL": 25.0, "Is.SPA_ancALL": False}))
    assert bool(m.iloc[0]) and "uncond" in which.iloc[0]


def test_conditioned_fallback_detected_without_an_isspa_c_column():
    # FELIX emits no Is.SPA_c; equality of p.value_c and p.value.NA_c is the flag
    m, which = fpr.spa_fallback(_frame(**{"Tstat_c_anc3": 30.0,
                                          "p.value_c_anc3": 1e-20,
                                          "p.value.NA_c_anc3": 1e-20}))
    assert bool(m.iloc[0]) and "cond" in which.iloc[0]


def test_low_score_is_not_a_fallback():
    m, _ = fpr.spa_fallback(_frame(**{"Tstat_ancALL": 0.5, "Is.SPA_ancALL": False}))
    assert not bool(m.iloc[0])


def test_case_carriers_sums_het_and_hom():
    assert int(fpr.case_carriers(_frame()).iloc[0]) == 12


# ---------------------------------------------------------------- bundle I/O

def _synthetic_bundle(tmp_path):
    """A minimal felix_run.summary.tar.gz with the members the module reads."""
    tmp_path = Path(tmp_path); tmp_path.mkdir(parents=True, exist_ok=True)
    tar = tmp_path / "felix_run.summary.tar.gz"
    cal = ("arm\tphenotype\tjoint_p_col\tlambda_gc_cct\tn_variants\tn_joint_tested\t"
           "n_gw_sig\tn_samples\tn_cases\n"
           "{arm}\tPH\tP_cct_admixed\t0.95\t1000\t1000\t{nsig}\t1000\t100\n")
    hits = _frame(**{"Tstat_ancALL": 25.0, "Is.SPA_ancALL": False})
    hits["P_cct_admixed_c"] = 1e-12
    hits["P_cct_admixed"] = 1e-12
    hits["CHR"], hits["POS"] = "chr1", 100
    with tarfile.open(tar, "w:gz") as t:
        def add(name, data: bytes):
            info = tarfile.TarInfo(name); info.size = len(data)
            t.addfile(info, io.BytesIO(data))
        for arm, nsig in (("main", 0), ("ru", 1), ("ru_biallelic", 0)):
            add(f"./{arm}/calibration_summary.tsv",
                cal.format(arm=arm, nsig=nsig).encode())
            add(f"./{arm}/qc/PH/top_hits.tsv",
                hits.to_csv(sep="\t", index=False).encode())
        add("./repeat_qc/chr1.loci.summary.json", json.dumps(
            {"n_samples": 1000, "hap_missing_frac": 0.0,
             "hap_possible_duplicate_frac": 0.0, "hap_ancestry_inconsistent_frac": 0.0,
             "counts": {"loci_written": 42, "records_unassigned": 1}}).encode())
    return tar


def test_load_bundle_and_reports(tmp_path):
    b = fpr.load_bundle(_synthetic_bundle(tmp_path), cache_dir=str(tmp_path))
    cal = fpr.calibration_report(b)
    assert set(cal.arm) == set(fpr.ARMS)
    assert (cal.status == "ok").all()
    # N_eff = 4*N*phi*(1-phi) with 100/1000 cases -> 360
    assert int(cal.N_eff.iloc[0]) == 360

    hits = fpr.hit_report(b)
    assert len(hits) == 3                       # one hit in each arm's top_hits
    assert (hits.verdict == "ARTIFACT").all()   # all are SPA fallbacks

    rq = fpr.repeat_qc_report(b)
    assert len(rq) == 1 and int(rq.loci_written.iloc[0]) == 42


def test_coverage_note_fires_when_hits_exceed_the_file(tmp_path):
    b = fpr.load_bundle(_synthetic_bundle(tmp_path), cache_dir=str(tmp_path))
    b["arms"]["ru"]["calibration"].loc[0, "n_gw_sig"] = 99   # file holds only 1
    notes = fpr.coverage_note(b)
    assert any("NOT screened" in n for n in notes)


def test_fetch_passes_through_local_paths(tmp_path):
    p = tmp_path / "x.txt"; p.write_text("hi")
    assert fpr.fetch(str(p)) == Path(str(p))


# ---------------------------------------------------------------- sweep tests

def test_cluster_positions_matches_arms_with_different_marker_ids():
    # the SNV arm and the repeat arms name the same site differently
    h = pd.DataFrame({"CHR": ["chr9", "chr9", "chr9", "chr1"],
                      "POS": [8991058, 8991058, 8991200, 500],
                      "MarkerID": ["chr9_8991058_T_TA", "9-8991058-8991069-A",
                                   "9-8991200-8991210-T", "1-500-510-A"]})
    lab = fpr.cluster_positions(h, window=1000)
    assert lab.iloc[0] == lab.iloc[1] == lab.iloc[2]   # within 1 kb -> one locus
    assert lab.iloc[3] != lab.iloc[0]                  # different chromosome


def test_cluster_positions_splits_beyond_the_window():
    h = pd.DataFrame({"CHR": ["chr1", "chr1"], "POS": [1000, 9000],
                      "MarkerID": ["a", "b"]})
    assert fpr.cluster_positions(h, window=1000).nunique() == 2


def _hits(rows):
    return pd.DataFrame(rows, columns=["run", "phenotype", "arm", "MarkerID",
                                       "CHR", "POS", "p", "case_carriers", "verdict"])


def test_corroborated_finds_one_locus_across_three_arms():
    h = _hits([("r", "PH", "main", "chr9_8991058_T_TA", "chr9", 8991058, 1e-9, 38, "keep"),
               ("r", "PH", "ru", "9-8991058-8991069-A", "chr9", 8991058, 1.4e-9, 38, "keep"),
               ("r", "PH", "ru_biallelic", "9-8991058-8991069-A", "chr9", 8991058,
                1.4e-9, 38, "keep")])
    c = fpr.corroborated(h)
    assert len(c) == 1 and int(c.n_arms.iloc[0]) == 3


def test_corroborated_ignores_artifacts():
    h = _hits([("r", "PH", "main", "m1", "chr1", 100, 1e-9, 30, "ARTIFACT"),
               ("r", "PH", "ru", "m2", "chr1", 100, 1e-9, 0, "ARTIFACT")])
    assert len(fpr.corroborated(h)) == 0


def test_recurrence_counts_distinct_phenotypes_not_rows():
    # same locus, same phenotype, three arms -> NOT recurrence
    h = _hits([("r", "A", "main", "m", "chr1", 100, 1e-9, 30, "keep"),
               ("r", "A", "ru", "m", "chr1", 100, 1e-9, 30, "keep"),
               ("r", "A", "ru_biallelic", "m", "chr1", 100, 1e-9, 30, "keep")])
    assert len(fpr.recurrence(h, min_phenotypes=2)) == 0
    h2 = _hits([("r1", "A", "ru", "m", "chr1", 100, 1e-9, 30, "keep"),
                ("r2", "B", "ru", "m", "chr1", 100, 1e-9, 30, "keep")])
    r = fpr.recurrence(h2, min_phenotypes=2)
    assert len(r) == 1 and int(r.n_phenotypes.iloc[0]) == 2


def test_sweep_over_local_bundles(tmp_path):
    from test_felix_phenotype_report import _synthetic_bundle  # noqa
    tars = {f"run{i}": str(_synthetic_bundle(tmp_path / f"b{i}")) for i in range(2)}
    tbl = pd.DataFrame([{"entity:x_id": k, "run_summary_tar": v} for k, v in tars.items()])
    cal, hits, fails = fpr.sweep(tbl, "entity:x_id", store=str(tmp_path / "store"),
                                 workers=2, progress=False)
    assert len(fails) == 0
    assert set(cal.run) == {"run0", "run1"}
    assert len(cal) == 6                       # 2 runs x 3 arms
    assert len(hits) == 6 and (hits.verdict == "ARTIFACT").all()


def test_harvest_is_cached(tmp_path):
    from test_felix_phenotype_report import _synthetic_bundle  # noqa
    tar = str(_synthetic_bundle(tmp_path))
    d1 = fpr.harvest(tar, "r", store=str(tmp_path / "s"))
    marker = d1 / "main" / "calibration_summary.tsv"
    marker.write_text("sentinel")             # second harvest must not overwrite
    d2 = fpr.harvest(tar, "r", store=str(tmp_path / "s"))
    assert d1 == d2 and marker.read_text() == "sentinel"


def test_sweep_notebook_always_refreshes_the_module():
    """Regression: felix_08 must not skip the bucket copy when the file exists.

    felix_07 writes felix_phenotype_report.py into the working directory via
    %%writefile, so a "fetch only if absent" guard silently imports whichever
    copy happens to be there — which produced
    `AttributeError: module has no attribute 'sweep'` against a pre-a67ebf5 copy.
    """
    nb_path = (Path(__file__).resolve().parents[2]
               / "notebooks" / "terra" / "felix_08_hit_sweep.ipynb")
    assert nb_path.exists(), f"notebook not found at {nb_path}"
    nb = json.loads(nb_path.read_text())
    setup = [ "".join(c["source"]) for c in nb["cells"]
              if c["cell_type"] == "code" and "MODULE = Path(" in "".join(c["source"])]
    assert len(setup) == 1, f"expected one bootstrap cell, found {len(setup)}"
    cell = setup[0]
    assert "if not MODULE.exists()" not in cell, \
        "bootstrap is conditional again — it must always re-copy from the bucket"
    assert "gsutil" in cell and "cp" in cell, "bootstrap no longer copies the module"
    # and it must verify what it imported rather than failing later on an attribute
    assert "sweep" in cell and "REQUIRED" in cell, \
        "bootstrap must assert the required functions are present"
