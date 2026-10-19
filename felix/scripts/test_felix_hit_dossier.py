"""Tests for felix_hit_dossier: the locus list, the derived checks, and the context annotation.

Run with the repo venv: ``venv/bin/python -m pytest felix/scripts/test_*.py -q``.

Nothing here touches the network -- ``locus_context`` takes an injectable ``fetch`` -- and the one test that reads the
real sweep table skips when the manuscript repo is not checked out beside this one.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import felix_hit_dossier as fhd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT_DATA = REPO_ROOT.parent / "aou-lr-phase-2-manuscript" / "data"


def hit(run, arm, chrom, pos, p, cc=50, verdict="keep", marker=None, **kw):
    return dict(run=run, arm=arm, CHR=f"chr{chrom}", POS=pos, p=p, case_carriers=cc, verdict=verdict,
                MarkerID=marker or f"chr{chrom}_{pos}_A_T", **kw)


# ------------------------------------------------------------------------------------------- the locus list


def test_retained_drops_excluded_and_thin_hits():
    df = pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1000, 1e-9),
        hit("felix.prod.A_1", "main", 1, 2000, 1e-9, verdict="EXCLUDE"),
        hit("felix.prod.A_1", "main", 1, 3000, 1e-9, cc=9),
        hit("felix.prod.A_1", "main", 23, 4000, 1e-9),  # chrX is not reviewed
    ])
    out = fhd.retained(df)
    assert list(out.POS) == [1000]


def test_retained_requires_the_screen_columns():
    with pytest.raises(KeyError, match="verdict"):
        fhd.retained(pd.DataFrame([{"run": "r", "arm": "main", "CHR": "chr1", "POS": 1, "p": 1e-9,
                                    "case_carriers": 50}]))


def test_loci_merge_within_the_window_and_split_beyond_it():
    df = pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
        hit("felix.prod.A_1", "main", 1, 1_400_000, 1e-10),  # 400 kb away: same locus, and the stronger lead
        hit("felix.prod.A_1", "main", 1, 2_000_000, 1e-9),  # 600 kb from the last: a second locus
    ])
    out = fhd.dossier_loci(df)
    assert len(out) == 2
    assert list(out.pos) == [1_400_000, 2_000_000]
    assert out.iloc[0].n_markers == 2
    assert out.iloc[0].p == 1e-10


def test_loci_separate_phenotypes_and_chromosomes_and_collect_arms():
    df = pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
        hit("felix.prod.A_1", "ru", 1, 1_000_100, 2e-9),
        hit("felix.prod.B_2", "main", 1, 1_000_000, 1e-9),  # same place, other phenotype
        hit("felix.prod.A_1", "main", 2, 1_000_000, 1e-9),  # same place, other chromosome
    ])
    out = fhd.dossier_loci(df)
    assert len(out) == 3
    a1 = out[(out.phenotype == "A_1") & (out.chrom == 1)].iloc[0]
    assert a1.arms == "main,ru" and a1.n_arms == 2 and a1.lead_arm == "main"


def test_locus_window_brackets_the_lead():
    df = pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9)])
    L = fhd.dossier_loci(df).iloc[0]
    assert L.window_start == 1_000_000 - fhd.REGION_HALF
    assert L.window_end == 1_000_000 + fhd.REGION_HALF
    assert L.locus_id == "A_1__chr1_1000000"


def test_window_start_never_goes_below_one():
    L = fhd.dossier_loci(pd.DataFrame([hit("felix.prod.A_1", "main", 1, 500, 1e-9)])).iloc[0]
    assert L.window_start == 1


@pytest.mark.parametrize("run,expected", [
    ("felix.prod.BI_174_2", "BI_174.2"), ("felix.prod.BI_174", "BI_174"), ("CV_413_11", "CV_413.11"),
])
def test_phecode_matches_the_figure_loader(run, expected):
    assert fhd.phecode(run) == expected


@pytest.mark.skipif(not (MANUSCRIPT_DATA / "felix_sweep_all_hits.tsv").is_file(),
                    reason="the manuscript repo's data/ is not available")
def test_locus_list_agrees_with_the_figure_3_loader():
    """The plate must review exactly the hits Figure 3 counts, so the two independent merges must give one answer."""
    sys.path.insert(0, str(REPO_ROOT))
    from manuscript_figures import fig3_real_data as frd

    hits = pd.read_csv(MANUSCRIPT_DATA / "felix_sweep_all_hits.tsv", sep="\t")
    mine = sorted((r.phenotype, r.chrom, r.pos) for r in fhd.dossier_loci(hits).itertuples())
    theirs = sorted((r["phenotype"], r["chrom"], r["pos"]) for r in frd.load().loci())
    assert mine == theirs


# ------------------------------------------------------------------------------------------- Hardy-Weinberg


def test_hardy_weinberg_on_a_locus_in_equilibrium():
    # 10,000 samples, allele frequency 0.2: 400 hom-alt, 3,200 het expected.
    out = fhd.hardy_weinberg(n_het=3200, n_hom=400, n_hap=20000)
    assert out["n_samples"] == 10000
    assert out["af"] == pytest.approx(0.2)
    assert out["F"] == pytest.approx(0.0, abs=1e-9)
    assert out["p"] == pytest.approx(1.0, abs=1e-6)


def test_hardy_weinberg_catches_the_fixed_heterozygous_record():
    """Every sample heterozygous: allele frequency 0.5 and F close to -1, the collapsed-paralogue signature."""
    out = fhd.hardy_weinberg(n_het=9900, n_hom=0, n_hap=20000)
    assert out["af"] == pytest.approx(0.495, abs=1e-3)
    assert out["F"] < -0.95
    assert out["het_frac"] > 0.98


def test_hardy_weinberg_returns_nan_when_it_cannot_run():
    for args in [(math.nan, 0, 1000), (10, 5, 0), (0, 0, 1000)]:  # no counts, no samples, monomorphic
        assert math.isnan(fhd.hardy_weinberg(*args)["F"])


def test_chi2_sf_matches_scipy():
    scipy_stats = pytest.importorskip("scipy.stats")
    for x, df in [(0.5, 1), (3.84, 1), (11.07, 5), (60.0, 5), (1.0, 2), (200.0, 5)]:
        assert fhd.chi2_sf(x, df) == pytest.approx(float(scipy_stats.chi2.sf(x, df)), rel=1e-9, abs=1e-300)


def test_chi2_isf_inverts_chi2_sf():
    for df in (1, 5):
        for p in (0.5, 1e-3, 1e-20):
            assert fhd.chi2_sf(fhd._chi2_isf(p, df), df) == pytest.approx(p, rel=1e-6)


# --------------------------------------------------------------------------------------------- lead records


def lead(betas, ses, ps, spa=False, p_het=None, p_hom=0.5, **kw):
    """A lead row with the per-ancestry block filled from three equal-length lists (anc1..5), plus ancALL."""
    rec = {"P_het_admixed_c": p_het, "P_hom_admixed_c": p_hom, "P_cct_admixed_c": p_het}
    for i, s in enumerate(fhd.STRATA[:-1]):
        rec[f"BETA_c_{s}"], rec[f"SE_c_{s}"], rec[f"p.value_c_{s}"] = betas[i], ses[i], ps[i]
        rec[f"Is.SPA_{s}"] = spa
        rec[f"Tstat_c_{s}"] = betas[i] / ses[i] ** 2 if ses[i] else 0.0
        rec[f"var_c_{s}"] = 1.0 / ses[i] ** 2 if ses[i] else 0.0
        rec[f"MissingRate_{s}"] = 0.0
    rec.update(kw)
    return pd.Series(rec)


def test_ancestry_block_has_one_row_per_stratum_with_flare_labels():
    blk = fhd.ancestry_block(lead([0.5] * 5, [0.1] * 5, [1e-6] * 5))
    assert list(blk.ancestry) == list(fhd.FLARE_ORDER) + ["ALL"]
    assert blk.iloc[0].z == pytest.approx(5.0)


def test_ancestry_block_sums_carriers_from_het_and_hom():
    blk = fhd.ancestry_block(pd.Series({"N_case_het_anc1": 7, "N_case_hom_anc1": 3}))
    assert blk.iloc[0].case_carriers == 10


def test_joint_coherence_is_about_one_for_a_healthy_row():
    """Five strata at z = 3 sum to chi2 45; a P_het drawn from that same chi2 must give ratio 1."""
    z2 = 5 * 9.0
    p_het = fhd.chi2_sf(z2, 5)
    co = fhd.joint_coherence(lead([0.3] * 5, [0.1] * 5, [2.7e-3] * 5, p_het=p_het))
    assert co["z2_sum"] == pytest.approx(z2)
    assert co["ratio"] == pytest.approx(1.0, rel=1e-3)


def test_joint_coherence_catches_a_p_the_strata_do_not_support():
    """The CV_446.2 signature: per-ancestry z of about 2 but a joint p of 1e-25."""
    co = fhd.joint_coherence(lead([0.2] * 5, [0.1] * 5, [0.045] * 5, p_het=1e-25))
    assert co["ratio"] > fhd.COHERENCE_MAX


def test_effect_consistency_separates_one_effect_from_several():
    same = fhd.effect_consistency(lead([0.5] * 5, [0.1] * 5, [1e-6] * 5))
    assert same["Q"] == pytest.approx(0.0, abs=1e-9) and same["same_sign"] == 1.0
    split = fhd.effect_consistency(lead([0.5, -0.5, 0.5, -0.5, 0.5], [0.1] * 5, [1e-6] * 5))
    assert split["p"] < fhd.HET_Q_P and split["same_sign"] < 1.0


def test_effect_consistency_needs_two_strata():
    out = fhd.effect_consistency(pd.Series({"BETA_c_anc1": 0.5, "SE_c_anc1": 0.1}))
    assert out["n_strata"] == 1 and out["df"] == 0 and math.isnan(out["p"])


def test_score_coherence_ignores_saddlepoint_corrected_strata():
    out = fhd.score_p_coherence(lead([0.5] * 5, [0.1] * 5, [1e-30] * 5, spa=True))
    assert out["n_plain"] == 0 and out["n_spa"] == 5 and math.isnan(out["max_log10_gap"])


def test_score_coherence_flags_an_uncorrected_mismatch():
    # z = 5 implies p about 6e-7; the row claims 1e-30.
    out = fhd.score_p_coherence(lead([0.5] * 5, [0.1] * 5, [1e-30] * 5, spa=False))
    assert out["n_plain"] == 5 and out["max_log10_gap"] > 20


# ------------------------------------------------------------------------------------------ regional support


def region(positions, ps, p_col=fhd.P_COL):
    return pd.DataFrame({"POS": positions, p_col: ps})


def test_regional_support_counts_neighbours_near_the_lead():
    out = fhd.regional_support(region([1000, 1500, 90_000, 400_000], [1e-9, 1e-8, 1e-7, 1e-12]), 1000, 1e-9)
    assert out["available"] and out["n_in_window"] == 3  # 400 kb is outside LD_PARTNER_WINDOW
    assert out["n_support"] == 2 and out["best_neighbour_p"] == 1e-8 and out["best_neighbour_bp"] == 500


def test_regional_support_requires_a_neighbour_near_the_lead_in_p():
    out = fhd.regional_support(region([1000, 1500], [1e-9, 0.4]), 1000, 1e-9)
    assert out["available"] and out["n_support"] == 0


def test_a_genome_wide_significant_neighbour_counts_however_strong_the_lead():
    """At p 7e-17 the relative rule would demand neighbours below 7e-14 and call this a lone peak."""
    out = fhd.regional_support(region([1000, 1500, 2000], [7.3e-17, 1e-12, 1e-10]), 1000, 7.3e-17)
    assert out["n_support"] == 2 and out["n_gw_sig"] == 2
    assert out["best_neighbour_p"] == 1e-12


def test_regional_support_is_unavailable_without_an_extract():
    assert fhd.regional_support(None, 1000, 1e-9)["available"] is False
    assert fhd.regional_support(pd.DataFrame(), 1000, 1e-9)["available"] is False


def test_an_extract_for_another_window_is_not_support():
    """A file cut for a different locus must read as 'not available', never as 'no neighbours'."""
    out = fhd.regional_support(region([5_000_000, 5_000_500], [1e-9, 1e-9]), 1000, 1e-9)
    assert out["available"] is False and out["n_markers"] == 2


# --------------------------------------------------------------------------------------------------- context


def fake_fetch(records):
    def fetch(track, chrom, start, end):
        return [r for r in records.get(track, [])
                if r.get("chromStart", r.get("txStart", 0)) < end and r.get("chromEnd", r.get("txEnd", 0)) > start]
    return fetch


def test_locus_context_reads_band_segdup_repeat_and_genes():
    fetch = fake_fetch({
        "cytoBand": [{"chrom": "chr9", "chromStart": 0, "chromEnd": 10_000_000, "name": "p24.1", "gieStain": "gneg"}],
        "genomicSuperDups": [{"chromStart": 999_000, "chromEnd": 1_001_000}],
        "simpleRepeat": [{"chromStart": 999_900, "chromEnd": 1_000_100, "period": 12, "copyNum": 8.0, "score": 100}],
        fhd.GENE_TRACK: [
            {"name2": "AAA", "txStart": 900_000, "txEnd": 1_100_000, "strand": "+", "chrom": "chr9"},
            {"name2": "BBB", "txStart": 1_200_000, "txEnd": 1_210_000, "strand": "-", "chrom": "chr9"},
        ],
    })
    ctx, genes = fhd.locus_context(9, 1_000_000, fetch=fetch)
    assert ctx["cytoband"] == "9p24.1" and ctx["pericentromeric"] is False
    assert ctx["in_segdup"] is True and ctx["segdup_frac"] > 0
    assert ctx["in_tandem_repeat"] is True and ctx["repeat_period"] == 12
    assert ctx["gene"] == "AAA" and ctx["nearest_gene"] == "AAA" and ctx["nearest_gene_bp"] == 0
    assert [g["gene"] for g in genes] == ["AAA", "BBB"]


def test_locus_context_reports_the_nearest_gene_when_the_lead_is_intergenic():
    fetch = fake_fetch({
        "cytoBand": [{"chrom": "chr1", "chromStart": 0, "chromEnd": 10_000_000, "name": "p36.13", "gieStain": "acen"}],
        fhd.GENE_TRACK: [{"name2": "CCC", "txStart": 1_050_000, "txEnd": 1_060_000, "strand": "+", "chrom": "chr1"}],
    })
    ctx, _ = fhd.locus_context(1, 1_000_000, fetch=fetch)
    assert ctx["gene"] == "" and ctx["nearest_gene"] == "CCC" and ctx["nearest_gene_bp"] == 50_000
    assert ctx["pericentromeric"] is True


def test_locus_context_merges_transcripts_of_one_gene():
    fetch = fake_fetch({fhd.GENE_TRACK: [
        {"name2": "DDD", "txStart": 990_000, "txEnd": 1_000_500, "strand": "+", "chrom": "chr1"},
        {"name2": "DDD", "txStart": 995_000, "txEnd": 1_020_000, "strand": "+", "chrom": "chr1"},
    ]})
    _, genes = fhd.locus_context(1, 1_000_000, fetch=fetch)
    assert len(genes) == 1 and genes[0]["start"] == 990_000 and genes[0]["end"] == 1_020_000


def test_annotate_context_survives_a_failed_lookup():
    def boom(*a, **k):
        raise TimeoutError("UCSC did not answer")

    loci = fhd.dossier_loci(pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9)]))
    ctx, genes, tracks = fhd.annotate_context(loci, fetch=boom, pause=0)
    assert len(ctx) == 1 and ctx.iloc[0].error.startswith("TimeoutError")
    assert genes.empty and tracks.empty


def test_locus_tracks_keeps_long_period_repeats_and_segdups():
    fetch = fake_fetch({
        "genomicSuperDups": [{"chromStart": 990_000, "chromEnd": 1_010_000}],
        "simpleRepeat": [
            {"chromStart": 999_000, "chromEnd": 1_001_000, "period": 28, "copyNum": 24.1, "score": 500},
            {"chromStart": 999_500, "chromEnd": 999_520, "period": 2, "copyNum": 10.0, "score": 50},
        ],
    })
    got = fhd.locus_tracks(9, 1_000_000, fetch=fetch)
    assert [t["kind"] for t in got] == ["repeat", "segdup"]
    assert got[0]["period"] == 28  # the 2 bp microsatellite is dropped


def test_locus_tracks_caps_each_kind():
    many = [{"chromStart": i * 100, "chromEnd": i * 100 + 50, "period": 20, "copyNum": 2.0, "score": 1}
            for i in range(200)]
    got = fhd.locus_tracks(1, 10_000, fetch=fake_fetch({"simpleRepeat": many}))
    assert len(got) == fhd.TRACK_MAX_PER_KIND


# ----------------------------------------------------------------------------------------- chips and verdict


def _status(flags, name):
    return next(f["status"] for f in flags if f["check"] == name)


def test_checks_mark_absent_inputs_not_available_rather_than_passing():
    L = fhd.dossier_loci(pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9)])).iloc[0]
    flags = fhd.checks(L, None, None, None)
    for name in ("Joint-test coherence", "Regional support", "Genomic context"):
        assert _status(flags, name) == "n/a"
    assert fhd.verdict(flags) == "incomplete"


def test_checks_fail_a_fixed_heterozygous_record():
    df = pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-70, ctrl_het_frac=0.996, ctrl_F=-0.993)])
    flags = fhd.checks(fhd.dossier_loci(df).iloc[0], None, None, None)
    assert _status(flags, "Hardy-Weinberg") == "fail"
    assert fhd.verdict(flags) == "fail"


def test_checks_warn_on_a_heterozygote_deficit():
    df = pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9, ctrl_het_frac=0.05, ctrl_F=0.45)])
    flags = fhd.checks(fhd.dossier_loci(df).iloc[0], None, None, None)
    assert _status(flags, "Hardy-Weinberg") == "warn"


def test_verdict_ranks_fail_above_warn_above_incomplete():
    assert fhd.verdict([{"status": "fail"}, {"status": "warn"}, {"status": "n/a"}]) == "fail"
    assert fhd.verdict([{"status": "warn"}, {"status": "n/a"}]) == "warn"
    assert fhd.verdict([{"status": "pass"}, {"status": "n/a"}]) == "incomplete"
    assert fhd.verdict([{"status": "pass"}]) == "pass"


# ------------------------------------------------------------------------------------------------- lead_rows


def test_lead_rows_takes_the_strongest_marker_in_the_window_per_arm():
    loci = fhd.dossier_loci(pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
        hit("felix.prod.A_1", "ru", 1, 1_000_000, 2e-9),
    ]))
    top = {
        ("felix.prod.A_1", "main"): pd.DataFrame({"POS": [1_000_000, 1_000_050, 9_000_000],
                                                  fhd.P_COL: [1e-9, 1e-12, 1e-40], "MarkerID": list("abc")}),
        ("felix.prod.A_1", "ru"): pd.DataFrame({"POS": [1_000_000], fhd.P_COL: [2e-9], "MarkerID": ["d"]}),
    }
    out = fhd.lead_rows(loci, top)
    main = out[out.arm == "main"].iloc[0]
    assert main.found and main.MarkerID == "b"  # strongest inside the window, not the 9 Mb row outside it
    assert bool(out[out.arm == "ru"].iloc[0].found)


def test_lead_rows_reports_an_arm_with_no_harvested_table():
    loci = fhd.dossier_loci(pd.DataFrame([hit("felix.prod.A_1", "ru", 1, 1_000_000, 1e-9)]))
    out = fhd.lead_rows(loci, {})
    assert not out.iloc[0].found and "no top_hits" in out.iloc[0].note


def test_region_plan_gives_one_span_per_locus_per_run():
    df = pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
        hit("felix.prod.A_1", "main", 1, 9_000_000, 1e-9),  # a second chr1 locus, far away
        hit("felix.prod.A_1", "main", 5, 1_000_000, 1e-9),
        hit("felix.prod.B_2", "main", 1, 1_000_000, 1e-9),
    ])
    plan = fhd.region_plan(fhd.dossier_loci(df))
    assert set(plan) == {"felix.prod.A_1", "felix.prod.B_2"}
    assert len(plan["felix.prod.A_1"]) == 3 and len(plan["felix.prod.B_2"]) == 1
    assert plan["felix.prod.A_1"][0]["start"] == 1_000_000 - fhd.REGION_HALF


def test_region_plan_merges_windows_that_overlap():
    """Reachable only when REGION_HALF is at least half LOCUS_WINDOW, so the loci frame is built directly here."""
    loci = pd.DataFrame([
        dict(locus_id="A__chr1_1000000", run="felix.prod.A_1", phenotype="A_1", chrom=1, pos=1_000_000,
             window_start=500_000, window_end=1_500_000),
        dict(locus_id="A__chr1_1600000", run="felix.prod.A_1", phenotype="A_1", chrom=1, pos=1_600_000,
             window_start=1_100_000, window_end=2_100_000),
    ])
    spans = fhd.region_plan(loci)["felix.prod.A_1"]
    assert len(spans) == 1
    assert spans[0]["start"] == 500_000 and spans[0]["end"] == 2_100_000
    assert spans[0]["locus_id"] == "A__chr1_1000000,A__chr1_1600000"


def test_region_plan_keeps_distant_windows_apart():
    df = pd.DataFrame([
        hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
        hit("felix.prod.A_1", "main", 1, 9_000_000, 1e-9),
    ])
    assert len(fhd.region_plan(fhd.dossier_loci(df))["felix.prod.A_1"]) == 2


def test_compact_region_drops_columns_and_thins_the_background(tmp_path):
    n = 5_000
    raw = pd.DataFrame({
        "CHR": ["chr1"] * n, "POS": range(900_000, 900_000 + n), "MarkerID": [f"m{i}" for i in range(n)],
        fhd.P_COL: [1e-9] * 5 + [0.4] * (n - 5),
        "junk_column": ["x"] * n,
    })
    src = tmp_path / "raw.tsv"
    raw.to_csv(src, sep="\t", index=False)
    locus = pd.Series({"window_start": 900_000, "window_end": 900_000 + n})
    out = tmp_path / "c.tsv.gz"
    written = fhd.compact_region(src, out, locus)
    got = pd.read_csv(out, sep="\t")
    assert written == 5 + fhd.REGION_BUDGET
    assert "junk_column" not in got.columns and "MarkerID" in got.columns
    assert (got[fhd.P_COL] < fhd.REGION_KEEP_BELOW).sum() == 5  # every strong marker survives


def test_compact_region_clips_to_the_locus_window(tmp_path):
    raw = pd.DataFrame({"CHR": ["chr1"] * 4, "POS": [10, 20, 30, 40], fhd.P_COL: [0.1] * 4})
    src = tmp_path / "raw.tsv"
    raw.to_csv(src, sep="\t", index=False)
    out = tmp_path / "c.tsv.gz"
    fhd.compact_region(src, out, pd.Series({"window_start": 20, "window_end": 30}))
    assert list(pd.read_csv(out, sep="\t").POS) == [20, 30]


def test_compact_region_refuses_a_table_without_the_p_column(tmp_path):
    src = tmp_path / "raw.tsv"
    pd.DataFrame({"CHR": ["chr1"], "POS": [1]}).to_csv(src, sep="\t", index=False)
    with pytest.raises(KeyError, match=fhd.P_COL):
        fhd.compact_region(src, tmp_path / "c.tsv.gz")


def test_thin_region_keeps_every_strong_row():
    df = pd.DataFrame({"POS": range(1000), fhd.P_COL: [1e-9] * 10 + [0.5] * 990})
    out = fhd._thin_region(df, fhd.P_COL, keep_below=1e-3, budget=100)
    assert (pd.to_numeric(out[fhd.P_COL]) < 1e-3).sum() == 10 and len(out) == 110


# -------------------------------------------------------------------------------------------- lead sites


def test_lead_sites_carries_alleles_for_a_main_arm_lead():
    df = pd.DataFrame([hit("felix.prod.A_1", "main", 9, 8_991_058, 1e-9, marker="chr9_8991058_T_TA")])
    out = fhd.lead_sites(fhd.dossier_loci(df)).iloc[0]
    assert (out.chrom, out.pos, out.ref, out.alt) == (9, 8_991_058, "T", "TA")


def test_lead_sites_carries_the_run_so_a_per_run_input_can_be_joined():
    """pheno_cov is a per-run input: a sex-restricted phenotype has its own file, so the run must travel along."""
    df = pd.DataFrame([hit("felix.prod.CA_142_2", "main", 9, 8_991_058, 1e-9)])
    out = fhd.lead_sites(fhd.dossier_loci(df)).iloc[0]
    assert out.run == "felix.prod.CA_142_2" and out.phenotype == "CA_142.2"


def test_lead_sites_leaves_alleles_blank_for_a_repeat_arm_lead():
    df = pd.DataFrame([hit("felix.prod.A_1", "ru", 9, 8_991_058, 1e-9, marker="9-8991058-8991069-A")])
    out = fhd.lead_sites(fhd.dossier_loci(df)).iloc[0]
    assert out.ref == "" and out.alt == "" and out.lead_arm == "ru"


def test_lead_sites_has_one_row_per_locus():
    df = pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9),
                       hit("felix.prod.A_1", "main", 1, 1_000_050, 2e-9),
                       hit("felix.prod.A_1", "main", 5, 1_000_000, 1e-9)])
    loci = fhd.dossier_loci(df)
    assert len(fhd.lead_sites(loci)) == len(loci) == 2


@pytest.mark.parametrize("marker,expected", [
    ("chr9_8991058_T_TA", ("T", "TA")),
    ("chr9_8991058_T_TA;chr9_8991058_TA_T;chr9_8991058_T_A_TA_TAA", ("TA", "TAA")),
    ("chr1:19786242:C:G", ("C", "G")),
    ("9-8991058-8991069-A", None),
    ("", None), (None, None),
])
def test_split_alleles_reads_the_last_two_fields(marker, expected):
    assert fhd.split_alleles(marker) == expected


@pytest.mark.skipif(not (MANUSCRIPT_DATA / "felix_sweep_all_hits.tsv").is_file(),
                    reason="the manuscript repo's data/ is not available")
def test_split_alleles_agrees_with_the_figure_loader():
    """This module duplicates the figure loader's allele split because it must import without the manuscript repo."""
    sys.path.insert(0, str(REPO_ROOT))
    from manuscript_figures import fig3_real_data as frd

    hits = pd.read_csv(MANUSCRIPT_DATA / "felix_sweep_all_hits.tsv", sep="\t")
    markers = [m for m in hits.MarkerID.dropna().unique()]
    assert len(markers) > 100
    assert all(fhd.split_alleles(m) == frd.split_alleles(m) for m in markers)


# ------------------------------------------------------------------------------------------- call rate


def _locus():
    return fhd.dossier_loci(pd.DataFrame([hit("felix.prod.A_1", "main", 1, 1_000_000, 1e-9)])).iloc[0]


def test_call_rate_is_not_assessed_without_the_pulled_genotypes():
    flags = fhd.checks(_locus(), None, None, None, None)
    chip = next(f for f in flags if f["check"] == "Call rate")
    assert chip["status"] == "n/a" and "pre-phasing" in chip["detail"]


def test_call_rate_passes_when_cases_and_controls_agree():
    cr = {"pre_call_rate_case": 0.991, "pre_call_rate_control": 0.994,
          "imputed_frac_case": 0.009, "imputed_frac_control": 0.006}
    chip = next(f for f in fhd.checks(_locus(), None, None, None, cr) if f["check"] == "Call rate")
    assert chip["status"] == "pass" and "0.991" in chip["detail"]


def test_call_rate_warns_on_a_case_control_gap():
    cr = {"pre_call_rate_case": 0.90, "pre_call_rate_control": 0.995,
          "imputed_frac_case": 0.10, "imputed_frac_control": 0.005}
    chip = next(f for f in fhd.checks(_locus(), None, None, None, cr) if f["check"] == "Call rate")
    assert chip["status"] == "warn" and "differs by" in chip["detail"]


def test_call_rate_warns_when_only_the_imputed_share_differs():
    """The failure mode a call-rate comparison alone cannot see: both well called, but filled in unequally."""
    cr = {"pre_call_rate_case": 0.99, "pre_call_rate_control": 0.99,
          "imputed_frac_case": 0.08, "imputed_frac_control": 0.01}
    chip = next(f for f in fhd.checks(_locus(), None, None, None, cr) if f["check"] == "Call rate")
    assert chip["status"] == "warn" and "imputed by phasing" in chip["detail"]


def test_call_rate_reports_a_note_when_the_locus_was_not_pulled():
    cr = {"note": "no A_1 column in pheno_cov"}
    chip = next(f for f in fhd.checks(_locus(), None, None, None, cr) if f["check"] == "Call rate")
    assert chip["status"] == "n/a" and "pheno_cov" in chip["detail"]


# ------------------------------------------------------------------------------- the Terra notebook itself

NOTEBOOK = REPO_ROOT / "notebooks" / "terra" / "felix_09_hit_dossier.ipynb"


@pytest.mark.skipif(not NOTEBOOK.is_file(), reason="the notebook is not in this checkout")
def test_every_notebook_code_cell_compiles():
    """A generated cell that does not parse only fails when it is run in Terra, several expensive cells in."""
    import ast
    import json

    nb = json.loads(NOTEBOOK.read_text())
    bad = []
    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        body = "".join(cell["source"])
        if body.lstrip().startswith("%%"):
            body = body.split("\n", 1)[1] if "\n" in body else ""
        body = "\n".join(l for l in body.split("\n") if not l.lstrip().startswith("!"))
        if not body.strip():
            continue
        try:
            ast.parse(body)
        except SyntaxError as exc:
            bad.append(f"cell {i}: {exc.msg} at line {exc.lineno}")
    assert not bad, "; ".join(bad)


@pytest.mark.skipif(not NOTEBOOK.is_file(), reason="the notebook is not in this checkout")
def test_the_notebook_bootstraps_every_module_it_imports():
    """A module imported in a later cell but absent from MODULES fails on Terra with ModuleNotFoundError."""
    import json
    import re

    nb = json.loads(NOTEBOOK.read_text())
    source = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    imported = set(re.findall(r"^\s*import (felix_\w+|extract_region|pull_site_fields)", source, re.M))
    setup = next("".join(c["source"]) for c in nb["cells"]
                 if c["cell_type"] == "code" and "MODULES" in "".join(c["source"]))
    for mod in imported:
        assert f"{mod}.py" in setup, f"{mod} is imported but not copied by the bootstrap"
        assert f'"{mod}"' in setup, f"{mod} is imported but not in the REQUIRED check"
    assert imported, "no local module imports found: the pattern needs updating"
