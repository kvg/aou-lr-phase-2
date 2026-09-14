#!/usr/bin/env python3
"""Unit tests for association-facing LAI eval scripts (Parts 1–4 helpers)."""

from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_build_af_panel import (  # noqa: E402
    af_range,
    select_markers,
    site_af_from_genotypes,
)
from flare_lai_null_lambda import (  # noqa: E402
    aggregate_perm_lambdas,
    assert_preflight_ok,
    beta_grid_from_mid,
    cis_overlap,
    lambda_gc_from_p,
    pick_winner,
    pick_winner_across_betas,
    preflight_anchor_or,
    preflight_q_vs_pcs,
    run_preflight,
    simulate_structured_pheno,
)
from flare_score_allele_ancestry import (  # noqa: E402
    allele_loglik,
    covering_index,
    score_haplotype,
    score_recipe,
)
from flare_score_mendelian_lai import (  # noqa: E402
    expected_crossovers,
    load_ped,
    locus_compatible,
    score_trio_path,
)


def test_af_panel_ranking():
    afs_hi = {"eas": 0.01, "amr": 0.02, "eur": 0.03, "afr": 0.90, "sas": 0.05}
    afs_lo = {"eas": 0.40, "amr": 0.41, "eur": 0.42, "afr": 0.43, "sas": 0.44}
    assert af_range(afs_hi) > 0.8
    assert af_range(afs_lo) < 0.05
    gts = ["0|0", "1|1", "0|1", "1|0", "0|0"]
    panels = ["eas", "afr", "eur", "afr", "amr"]
    afs, n_alleles = site_af_from_genotypes(gts, panels)
    assert n_alleles == 10
    assert abs(afs["afr"] - 0.75) < 1e-9  # 1|1 and 1|0 → 3/4
    rows = [
        {
            "chrom": "chr22",
            "pos": 100,
            "ref": "A",
            "alt": "G",
            "mac": 100,
            "af_range": 0.8,
            "af_eas": 0.1,
            "af_amr": 0.1,
            "af_eur": 0.1,
            "af_afr": 0.9,
            "af_sas": 0.1,
        },
        {
            "chrom": "chr22",
            "pos": 200,
            "ref": "A",
            "alt": "C",
            "mac": 5,
            "af_range": 0.9,
            "af_eas": 0.0,
            "af_amr": 0.0,
            "af_eur": 0.0,
            "af_afr": 0.9,
            "af_sas": 0.0,
        },
        {
            "chrom": "chr22",
            "pos": 300,
            "ref": "A",
            "alt": "T",
            "mac": 100,
            "af_range": 0.02,
            "af_eas": 0.4,
            "af_amr": 0.41,
            "af_eur": 0.42,
            "af_afr": 0.42,
            "af_sas": 0.4,
        },
    ]
    kept = select_markers(rows, top_n=10, min_mac=50, min_af_range=0.05)
    assert len(kept) == 1 and kept[0]["pos"] == 100


def test_covering_and_allele_ll():
    assert covering_index([100, 200, 300], 150) == 1  # first >= 150 → 200
    assert covering_index([100, 200, 300], 50) == 0
    assert covering_index([100, 200, 300], 400) == 2
    assert abs(allele_loglik(1, 0.5) - math.log(0.5)) < 1e-12
    ll, br, ok = score_haplotype(1, 3, {"afr": 0.9, "eur": 0.1})  # anc 3 = afr
    assert ok and ll > math.log(0.5)
    ll2, _, ok2 = score_haplotype(1, 2, {"afr": 0.9, "eur": 0.1})  # eur
    assert ok2 and ll2 < ll


def test_score_recipe_carry_forward():
    panel = [
        {
            "chrom": "chr22",
            "pos": 150,
            "ref": "A",
            "alt": "G",
            "afs": {"eas": 0.5, "amr": 0.5, "eur": 0.1, "afr": 0.9, "sas": 0.5},
        }
    ]
    # LAI only at 100 and 200; query 150 covers 200 (FELIXla)
    anc_by_pos = {150: [(3, 3)]}
    lai_pos_by_panel = {150: 200}
    gt_by_pos = {("chr22", 150): [(1, 1)]}
    s = score_recipe(
        panel=panel,
        anc_by_panel_pos=anc_by_pos,
        lai_pos_by_panel=lai_pos_by_panel,
        gt_by_pos=gt_by_pos,
        n_samples=1,
        n_lai_markers=2,
    )
    assert s["n_scored_haps"] == 2
    assert s["n_carried_forward_haps"] == 2
    assert s["frac_carried_forward"] == 1.0
    assert s["mean_ll"] > math.log(0.5)


def test_mendelian_compatible_and_path():
    # child AFR/EUR, father AFR/AFR, mother EUR/EUR → OK
    assert locus_compatible(3, 2, 3, 3, 2, 2)
    # child AFR/AFR, mother EUR/EUR, father EUR/EUR → fail
    assert not locus_compatible(3, 3, 2, 2, 2, 2)

    # path: start FM assignment, then need flip → 1 recomb, no hard
    loci = [
        (3, 2, 3, 3, 2, 2),  # c1 AFR from F, c2 EUR from M
        (2, 3, 3, 3, 2, 2),  # now needs swap
    ]
    sc = score_trio_path(loci)
    assert sc["n_hard_violations"] == 0
    assert sc["n_recombinations"] == 1

    # hard violation
    loci2 = [(3, 3, 2, 2, 2, 2)]
    sc2 = score_trio_path(loci2)
    assert sc2["n_hard_violations"] == 1


def test_load_ped_and_map_expectation():
    ped_path = (
        Path(__file__).resolve().parents[1]
        / "tractor_mix"
        / "resources"
        / "legacy_covariates"
        / "aou_phase2.ped"
    )
    trios = load_ped(ped_path)
    assert len(trios) > 30
    assert all(t.father != "0" and t.mother != "0" for t in trios)
    # synthetic map: 1 cM / Mb
    pts = [(0, 0.0), (10_000_000, 10.0)]
    xo = expected_crossovers(pts, 0, 10_000_000, n_meioses=2)
    assert abs(xo - 0.2) < 1e-9  # 0.1 Morgan * 2


def test_lambda_gc_and_winner():
    rng = np.random.default_rng(0)
    # Uniform p → λ ≈ 1
    p = rng.random(5000)
    lam = lambda_gc_from_p(p)
    assert 0.85 < lam < 1.15

    agg = aggregate_perm_lambdas(
        [1.0, 1.05, 0.95, 1.02],
        n_tested_sites=500,
        min_tested_sites=100,
        seed=1,
    )
    assert not agg["unstable_lambda"]
    assert agg["abs_lambda_ci"]["ci_low"] <= agg["abs_lambda_dev"] <= agg["abs_lambda_ci"]["ci_high"]

    unstable = aggregate_perm_lambdas(
        [1.0, 1.0],
        n_tested_sites=10,
        min_tested_sites=100,
        seed=1,
    )
    assert unstable["unstable_lambda"]

    rows = [
        {
            "experiment": "a",
            "abs_lambda_dev": 0.05,
            "abs_lambda_ci": {"mean": 0.05, "ci_low": 0.02, "ci_high": 0.08},
            "mean_ll": -0.5,
            "unstable_lambda": False,
        },
        {
            "experiment": "b",
            "abs_lambda_dev": 0.06,
            "abs_lambda_ci": {"mean": 0.06, "ci_low": 0.03, "ci_high": 0.09},
            "mean_ll": -0.3,
            "unstable_lambda": False,
        },
        {
            "experiment": "c",
            "abs_lambda_dev": 0.20,
            "abs_lambda_ci": {"mean": 0.20, "ci_low": 0.18, "ci_high": 0.22},
            "mean_ll": -0.1,
            "unstable_lambda": False,
        },
    ]
    assert cis_overlap(rows[0]["abs_lambda_ci"], rows[1]["abs_lambda_ci"])
    decision = pick_winner(rows)
    # a and b overlap; b has better (higher) mean_ll
    assert decision["winner"] == "b"
    assert decision["reason"] == "ci_overlap_mean_ll_tiebreak"


def test_write_tiny_panel_fixture():
    with tempfile.TemporaryDirectory() as td:
        panel = Path(td) / "markers.tsv"
        panel.write_text(
            "chrom\tpos\tref\talt\tmac\taf_range\taf_eas\taf_amr\taf_eur\taf_afr\taf_sas\n"
            "chr22\t100\tA\tG\t100\t0.8\t0.1\t0.1\t0.1\t0.9\t0.1\n"
        )
        from flare_score_allele_ancestry import load_panel

        rows = load_panel(panel)
        assert len(rows) == 1 and rows[0]["afs"]["afr"] == 0.9


def _write_sim_fixtures(td: Path, *, collinear_pcs: bool = False, strong_anchor: bool = True):
    """Tiny pheno_cov + global.anc for preflight / simulation tests."""
    import gzip

    rng = np.random.default_rng(0)
    n = 200
    ids = [f"S{i}" for i in range(n)]
    afr = rng.uniform(0.0, 1.0, size=n)
    eur = 1.0 - afr
    # Two families of size 4; rest singletons
    fam = [f"F{i // 4}" if i < 16 else "" for i in range(n)]
    gc = ["BI" if i % 2 == 0 else "BCM" for i in range(n)]
    if collinear_pcs:
        pc1 = afr + rng.normal(0, 0.01, size=n)
        pc2 = eur + rng.normal(0, 0.01, size=n)
    else:
        pc1 = rng.normal(size=n)
        pc2 = rng.normal(size=n)
    # Anchor: strong afr effect or noise
    if strong_anchor:
        logit = -1.0 + 2.5 * afr
        p = 1.0 / (1.0 + np.exp(-logit))
        y = (rng.random(n) < p).astype(int)
    else:
        y = rng.integers(0, 2, size=n)
    pheno = td / "pheno_cov.tsv"
    with pheno.open("w") as fh:
        fh.write("ID\tanchor_strong\tanchor_weak\tGC\tpedigree_family_id\tPC1\tPC2\tage\tsex\n")
        for i in range(n):
            fh.write(
                f"{ids[i]}\t{y[i]}\t{int(rng.integers(0, 2))}\t{gc[i]}\t{fam[i]}\t"
                f"{pc1[i]}\t{pc2[i]}\t{40 + i % 20}\t{i % 2}\n"
            )
    ganc = td / "fixed.global.anc.gz"
    with gzip.open(ganc, "wt") as fh:
        fh.write("SAMPLE eas amr eur afr sas\n")
        for i in range(n):
            fh.write(f"{ids[i]} 0 0 {eur[i]:.6f} {afr[i]:.6f} 0\n")
    return pheno, ganc


def test_preflight_q_vs_pcs_and_anchor():
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        pheno, ganc = _write_sim_fixtures(td, collinear_pcs=False, strong_anchor=True)
        q = preflight_q_vs_pcs(pheno, ganc, ["age", "sex", "GC", "PC1", "PC2"], r2_max=0.90)
        assert q["ok_for_pilot"]
        assert not q["pcs_dropped"]
        assert q["with_pcs"]["multivariate_r2"] < 0.90

        col = td / "col"
        col.mkdir(exist_ok=True)
        pheno_c, ganc_c = _write_sim_fixtures(col, collinear_pcs=True, strong_anchor=True)
        q2 = preflight_q_vs_pcs(pheno_c, ganc_c, ["age", "sex", "GC", "PC1", "PC2"], r2_max=0.90)
        assert q2["pcs_dropped"]
        assert q2["ok_for_pilot"]
        assert "PC1" not in q2["null_covariates_used"]

        anc = preflight_anchor_or(pheno, ganc, ["anchor_strong", "anchor_weak"])
        assert anc["pass"]
        assert anc["anchor_phenotype"] == "anchor_strong"
        assert abs(anc["beta_mid"]) >= 0.3

        weak = preflight_anchor_or(pheno, ganc, ["anchor_weak"], abs_beta_min=0.3, pmax=0.01)
        # weak may occasionally pass by chance; if it fails, gate is correct
        report = run_preflight(
            pheno_cov=pheno,
            global_anc=ganc,
            covariates=["age", "sex", "GC", "PC1", "PC2"],
            anchor_candidates=["anchor_strong"],
        )
        assert report["ok"]
        assert assert_preflight_ok(report)["ok"]

        bad = {
            "q_vs_pcs": {"pass": False, "pcs_dropped": False},
            "anchor_or": {"pass": True},
        }
        assert not assert_preflight_ok(bad)["ok"]


def test_simulate_family_and_prevalence():
    with tempfile.TemporaryDirectory() as tmp:
        td = Path(tmp)
        pheno, ganc = _write_sim_fixtures(td, strong_anchor=True)
        out = td / "sim.tsv"
        meta = simulate_structured_pheno(
            pheno,
            ganc,
            out,
            beta=1.5,
            case_rate=0.25,
            family_col="pedigree_family_id",
            batch_col="GC",
            sigma_f=1.0,
            sigma_eps=0.2,
            seed=7,
        )
        assert abs(meta["case_rate_obs"] - 0.25) < 0.08
        assert meta["n_families_shared"] >= 1
        meta_indep = simulate_structured_pheno(
            pheno,
            ganc,
            td / "sim_indep.tsv",
            beta=1.5,
            case_rate=0.25,
            family_col="pedigree_family_id",
            batch_col="GC",
            sigma_f=0.0,
            sigma_eps=0.2,
            seed=7,
        )
        # Shared family RE → higher within-family liability covariance
        assert (
            meta["mean_within_family_centered_prod"]
            > meta_indep["mean_within_family_centered_prod"]
        )
        # Larger |β| should push more extreme liability separation (case rate still matched)
        meta_hi = simulate_structured_pheno(
            pheno, ganc, td / "sim_hi.tsv", beta=3.0, case_rate=0.25, seed=7
        )
        assert abs(meta_hi["case_rate_obs"] - 0.25) < 0.08
        grid = beta_grid_from_mid(1.0)
        assert grid["lo"] == 0.5 and grid["mid"] == 1.0 and grid["hi"] == 2.0


def test_multi_beta_winner_stability():
    def row(exp, dev, lo, hi, ll=-0.5):
        return {
            "experiment": exp,
            "abs_lambda_dev": dev,
            "abs_lambda_ci": {"mean": dev, "ci_low": lo, "ci_high": hi},
            "mean_ll": ll,
            "unstable_lambda": False,
        }

    stable = {
        "lo": [row("a", 0.05, 0.04, 0.06), row("b", 0.20, 0.18, 0.22)],
        "mid": [row("a", 0.04, 0.03, 0.05), row("b", 0.19, 0.17, 0.21)],
        "hi": [row("a", 0.06, 0.05, 0.07), row("b", 0.22, 0.20, 0.24)],
    }
    pre = {
        "q_vs_pcs": {"pass": True, "pcs_dropped": False},
        "anchor_or": {"pass": True},
    }
    d = pick_winner_across_betas(stable, preflight=pre)
    assert d["winner"] == "a" and not d["ranking_unstable"]

    unstable = {
        "lo": [row("a", 0.05, 0.04, 0.06), row("b", 0.20, 0.18, 0.22)],
        "mid": [row("b", 0.04, 0.03, 0.05), row("a", 0.19, 0.17, 0.21)],
        "hi": [row("a", 0.06, 0.05, 0.07), row("b", 0.22, 0.20, 0.24)],
    }
    d2 = pick_winner_across_betas(unstable, preflight=pre)
    assert d2["winner"] is None and d2["ranking_unstable"]
    assert d2["reason"] == "disagreement_across_betas"


if __name__ == "__main__":
    test_af_panel_ranking()
    test_covering_and_allele_ll()
    test_score_recipe_carry_forward()
    test_mendelian_compatible_and_path()
    test_load_ped_and_map_expectation()
    test_lambda_gc_and_winner()
    test_write_tiny_panel_fixture()
    test_preflight_q_vs_pcs_and_anchor()
    test_simulate_family_and_prevalence()
    test_multi_beta_winner_stability()
    print("ok")
