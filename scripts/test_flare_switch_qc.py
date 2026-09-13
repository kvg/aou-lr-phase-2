#!/usr/bin/env python3
"""Unit tests for flare_switch_qc (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_switch_qc import (  # noqa: E402
    SampleTracker,
    SwitchEvent,
    assign_nearest_het,
    build_tracts,
    call_qc_filtered_tract_metrics,
    classify_gt,
    expected_mean_tract_mb,
    expected_switches_per_hap_per_mb,
    filter_switches_call_qc,
    implied_t_gen,
    is_glnexus_homref_sentinel,
    mark_flickers,
    maybe_record_switch,
    process_sample_site,
    rnc_has_I,
    summarize,
    switch_fails_call_qc,
    QcAccum,
    GtBin,
    parse_sample_fields,
    query_format_string,
    sample_format_tags,
)


def test_detects_switch_and_flicker():
    tr = SampleTracker()
    maybe_record_switch(tr, "S1", 1, "chr22", 100, "A", "G", 2, 40, 20, "..")
    maybe_record_switch(tr, "S1", 1, "chr22", 200, "C", "T", 3, 5, 2, "II")
    maybe_record_switch(tr, "S1", 1, "chr22", 300, "G", "A", 2, 30, 15, "..")
    assert len(tr.switches) == 2
    mark_flickers(tr.switches, max_span_bp=50_000)
    assert tr.switches[0].is_flicker
    assert tr.switches[1].is_flicker
    assert tr.switches[0].anc_from == 2
    assert tr.switches[0].anc_to == 3


def test_no_switch_when_anc_stable():
    tr = SampleTracker()
    maybe_record_switch(tr, "S1", 1, "chr22", 100, "A", "G", 2, 40, 20, "..")
    maybe_record_switch(tr, "S1", 1, "chr22", 200, "C", "T", 2, 40, 20, "..")
    assert tr.switches == []


def test_summarize_enrichment():
    events = [
        SwitchEvent(
            sample="S1",
            hap=1,
            chrom="chr22",
            pos=100,
            ref="A",
            alt="G",
            anc_from=2,
            anc_to=3,
            gq=5,
            dp=2,
            rnc="II",
            prev_pos=50,
            bp_gap=50,
            is_flicker=True,
            flicker_span_bp=50,
        )
    ]
    qc = QcAccum(
        switch_gq=[5],
        switch_dp=[2],
        bg_gq=[40, 40, 40],
        bg_dp=[20, 20, 20],
        switch_low_gq=1,
        switch_low_dp=1,
        switch_rnc_i=1,
        switch_n=1,
        bg_low_gq=0,
        bg_low_dp=0,
        bg_rnc_i=0,
        bg_n=3,
    )
    summary = summarize(
        qc,
        {"gq_threshold": 20, "dp_threshold": 5, "n_samples": 1, "sites_scanned": 10},
        events,
    )
    assert summary["enrichment_frac_low_gq"] is None or summary["enrichment_frac_low_gq"] > 1
    # bg low-gq frac is 0 → enrichment None
    assert summary["switch_site_calls"]["frac_low_gq"] == 1.0
    assert summary["switches_flicker"]["n"] == 1
    assert "het" in summary["by_gt"]
    assert summary["switch_site_calls"]["frac_sentinel"] == 0.0


def test_classify_gt_and_sentinel():
    assert classify_gt("0|0") == "hom_ref"
    assert classify_gt("0/0") == "hom_ref"
    assert classify_gt("0|1") == "het"
    assert classify_gt("1/0") == "het"
    assert classify_gt("1|1") == "hom_alt"
    assert classify_gt(".") == "missing"
    assert classify_gt("./.") == "missing"
    assert is_glnexus_homref_sentinel("hom_ref", 1, 0)
    assert is_glnexus_homref_sentinel("hom_ref", 0, 0)
    assert not is_glnexus_homref_sentinel("hom_ref", 30, 0)
    assert not is_glnexus_homref_sentinel("het", 1, 0)


def test_summarize_het_enrichment_not_confounded_by_homref():
    qc = QcAccum()
    qc.switch_n = 10
    qc.bg_n = 10
    qc.switch_low_gq = 8
    qc.bg_low_gq = 8
    qc.switch_by_gt["hom_ref"] = GtBin(n=8, low_gq=8, sentinel=8, gq=[1] * 8, dp=[0] * 8)
    qc.switch_by_gt["het"] = GtBin(n=2, low_gq=0, sentinel=0, gq=[40, 45], dp=[12, 15])
    qc.bg_by_gt["hom_ref"] = GtBin(n=8, low_gq=8, sentinel=8, gq=[1] * 8, dp=[0] * 8)
    qc.bg_by_gt["het"] = GtBin(n=2, low_gq=0, sentinel=0, gq=[40, 42], dp=[14, 16])
    summary = summarize(
        qc,
        {"gq_threshold": 20, "dp_threshold": 5, "n_samples": 1, "sites_scanned": 10},
        [],
    )
    assert summary["by_gt"]["het"]["enrichment_frac_low_gq"] is None
    assert summary["by_gt"]["hom_ref"]["switch"]["frac_sentinel"] == 1.0
    assert summary["by_gt"]["het"]["switch"]["median_gq"] == 42.5


def _site(tr, pos, an, gq, dp, gt, rnc=".."):
    return process_sample_site(
        tr, "S1", "chr22", pos, "A", "G", an, an, gq, dp, rnc, gt,
    )


def test_nearest_het_for_homref_switch():
    tr = SampleTracker()
    _site(tr, 100, 2, 40, 18, "0|1")
    ev1, ev2, gt_class, _ = _site(tr, 200, 3, 1, 0, "0|0")
    assert gt_class == "hom_ref"
    assert ev1 is not None
    _site(tr, 230, 3, 8, 3, "0|1")
    assign_nearest_het(ev1)
    assert ev1.prev_het_pos == 100
    assert ev1.prev_het_gq == 40
    assert ev1.prev_het_gap == 100
    assert ev1.next_het_pos == 230
    assert ev1.next_het_gq == 8
    assert ev1.next_het_gap == 30
    assert ev1.nearest_het_side == "next"
    assert ev1.nearest_het_gq == 8
    assert ev1.nearest_het_gap == 30
    assert ev1.nearest_het_rnc == ".."
    assert ev1.prev_het_rnc == ".."


def test_nearest_het_is_current_when_switch_is_het():
    tr = SampleTracker()
    _site(tr, 100, 2, 40, 18, "0|1")
    ev1, _, gt_class, _ = _site(tr, 200, 3, 5, 2, "0|1")
    _site(tr, 400, 3, 50, 20, "0|1")
    assert gt_class == "het"
    assert ev1 is not None
    assign_nearest_het(ev1)
    assert ev1.nearest_het_side == "current"
    assert ev1.nearest_het_gq == 5
    assert ev1.nearest_het_gap == 0
    assert ev1.prev_het_pos == 100
    assert ev1.next_het_pos == 400


def test_prev_het_does_not_cross_chromosomes():
    tr = SampleTracker()
    _site(tr, 100, 2, 40, 18, "0|1")
    process_sample_site(tr, "S1", "chr21", 50, "A", "G", 2, 2, 30, 12, "..", "0|1")
    ev1, _, _, _ = process_sample_site(
        tr, "S1", "chr21", 80, "C", "T", 3, 3, 1, 0, "..", "0|0",
    )
    assert ev1 is not None
    assign_nearest_het(ev1)
    assert ev1.prev_het_pos == 50
    assert ev1.prev_het_gq == 30


def test_summarize_nearest_het_enrichment():
    events = [
        SwitchEvent(
            sample="S1",
            hap=1,
            chrom="chr22",
            pos=200,
            ref="A",
            alt="G",
            anc_from=2,
            anc_to=3,
            gq=1,
            dp=0,
            rnc="..",
            prev_pos=100,
            bp_gap=100,
            gt="0|0",
            gt_class="hom_ref",
            nearest_het_gq=5,
            nearest_het_dp=2,
            nearest_het_gap=30,
            nearest_het_side="next",
        )
    ]
    qc = QcAccum()
    qc.bg_by_gt["het"] = GtBin(n=10, low_gq=1, low_dp=1, gq=[40] * 10, dp=[16] * 10)
    summary = summarize(
        qc,
        {"gq_threshold": 20, "dp_threshold": 5, "n_samples": 1, "sites_scanned": 10},
        events,
    )
    nh = summary["nearest_het"]
    assert nh["all"]["median_gq"] == 5
    assert nh["switch_at_hom_ref"]["n"] == 1
    assert nh["enrichment_frac_low_gq"] == 1.0 / 0.1
    assert nh["side_counts"]["next"] == 1
    assert "rnc" in summary


def test_homref_switch_records_incomplete_rnc():
    tr = SampleTracker()
    _site(tr, 100, 2, 40, 18, "0|1", rnc="..")
    ev1, _, _, _ = _site(tr, 200, 3, 1, 0, "0|0", rnc="II")
    assert ev1 is not None
    assert ev1.rnc == "II"
    assert rnc_has_I(ev1.rnc)
    _site(tr, 250, 3, 8, 3, "0|1", rnc="I.")
    assign_nearest_het(ev1)
    assert ev1.nearest_het_rnc == "I."
    assert ev1.prev_het_rnc == ".."


def test_summarize_rnc_enrichment():
    events = [
        SwitchEvent(
            sample="S1",
            hap=1,
            chrom="chr22",
            pos=200,
            ref="A",
            alt="G",
            anc_from=2,
            anc_to=3,
            gq=1,
            dp=0,
            rnc="II",
            prev_pos=100,
            bp_gap=100,
            gt="0|0",
            gt_class="hom_ref",
        )
    ]
    qc = QcAccum()
    qc.switch_n = 1
    qc.switch_rnc_i = 1
    qc.switch_rnc = ["II"]
    qc.bg_n = 10
    qc.bg_rnc_i = 1
    qc.bg_rnc = [".."] * 9 + ["II"]
    qc.bg_by_gt["het"] = GtBin(n=10, low_gq=1, low_dp=1, gq=[40] * 10, dp=[16] * 10)
    summary = summarize(
        qc,
        {"gq_threshold": 20, "dp_threshold": 5, "n_samples": 1, "sites_scanned": 10},
        events,
    )
    assert summary["rnc"]["switch"]["frac_I"] == 1.0
    assert summary["rnc"]["background"]["frac_I"] == 0.1
    assert summary["rnc"]["enrichment_frac_I"] == 10.0
    assert summary["rnc"]["switch_at_hom_ref"]["frac_I"] == 1.0
    assert rnc_has_I("II")
    assert rnc_has_I("I.")
    assert not rnc_has_I("..")


def test_expected_rate_from_t():
    assert expected_mean_tract_mb(100.0) == 1.0
    # Without props: factor=1 → rate = T/100
    assert expected_switches_per_hap_per_mb(100.0) == 1.0
    assert expected_mean_tract_mb(6.0) == 100.0 / 6.0
    assert abs(implied_t_gen(1.2) - 120.0) < 1e-9
    # Two-way 50/50: 1 - sum p^2 = 0.5
    rate = expected_switches_per_hap_per_mb(100.0, [0.5, 0.5])
    assert abs(rate - 0.5) < 1e-9
    from flare_switch_qc import het_admixture_factor, implied_T_given_props, span_ok_for_t

    assert abs(implied_T_given_props(0.5, [0.5, 0.5]) - 100.0) < 1e-9
    # mean tract at T=8 is 12.5 Mb; 3*12.5=37.5 → 30 Mb is False
    assert span_ok_for_t(30.0, 8.0) is False
    assert span_ok_for_t(40.0, 8.0) is True
    assert abs(het_admixture_factor([0.8, 0.2]) - (1 - 0.64 - 0.04)) < 1e-12


def test_build_tracts_includes_no_switch_haps_and_censoring():
    events = [
        SwitchEvent(
            sample="S1",
            hap=1,
            chrom="chr22",
            pos=200,
            ref="A",
            alt="G",
            anc_from=2,
            anc_to=3,
            gq=40,
            dp=12,
            rnc="..",
            prev_pos=100,
            bp_gap=100,
        ),
        SwitchEvent(
            sample="S1",
            hap=1,
            chrom="chr22",
            pos=250,
            ref="C",
            alt="T",
            anc_from=3,
            anc_to=2,
            gq=5,
            dp=2,
            rnc="..",
            prev_pos=200,
            bp_gap=50,
            is_flicker=True,
            flicker_span_bp=50,
        ),
    ]
    events[0].is_flicker = True
    events[0].flicker_span_bp = 50
    tracts = build_tracts(
        events,
        samples=["S1", "S2"],
        last_anc={("S1", 1): 2, ("S1", 2): 2, ("S2", 1): 3, ("S2", 2): 3},
        pos_min={"chr22": 100},
        pos_max={"chr22": 400},
    )
    s2 = [t for t in tracts if t.sample == "S2" and t.hap == 1][0]
    assert s2.length_bp == 300
    assert s2.left_censored and s2.right_censored
    assert s2.n_switches_on_hap == 0
    s1 = [t for t in tracts if t.sample == "S1" and t.hap == 1]
    assert [t.start for t in s1] == [100, 200, 250]
    assert s1[0].left_censored and not s1[0].right_censored
    assert not s1[1].left_censored and not s1[1].right_censored
    assert s1[1].is_flicker_interval
    assert s1[1].length_bp == 50
    assert s1[2].right_censored
    summary = summarize(
        QcAccum(),
        {
            "gq_threshold": 20,
            "dp_threshold": 5,
            "n_samples": 2,
            "sites_scanned": 10,
            "span_bp": 300,
        },
        events,
        tracts,
    )
    tr = summary["tracts"]
    assert tr["n_haps"] == 4
    assert tr["n_haps_with_switch"] == 1
    assert abs(tr["switches_per_hap_per_mb"] - (2 / 4 / (300 / 1e6))) < 1e-6
    assert tr["flicker_frac"] == 1.0


def test_an_only_query_omits_gq_dp():
    tags = sample_format_tags(
        {"AN1": True, "AN2": True, "GT": True, "GQ": True, "DP": True, "RNC": True},
        an_only=True,
    )
    assert tags == ["AN1", "AN2"]
    fmt = query_format_string(tags, an_only=True)
    assert "%GQ" not in fmt
    assert "%REF" not in fmt
    assert "%AN1" in fmt
    tags_full = sample_format_tags(
        {"AN1": True, "AN2": True, "GT": True, "GQ": False, "DP": False, "RNC": False}
    )
    assert tags_full == ["AN1", "AN2", "GT"]
    parts = ["chr22", "100", "A", "G", "3", "2", "0|1"]
    gq, dp, an1, an2, rnc, gt = parse_sample_fields(parts, 4, tags_full)
    assert gq is None and dp is None
    assert an1 == 3 and an2 == 2
    assert rnc == "." and gt == "0|1"


def test_split_sample_chunks():
    from flare_switch_qc import split_sample_chunks

    assert split_sample_chunks(["a", "b", "c", "d"], 1) == [["a", "b", "c", "d"]]
    assert split_sample_chunks(["a", "b", "c", "d"], 2) == [["a", "b"], ["c", "d"]]
    assert split_sample_chunks(["a", "b", "c"], 2) == [["a", "b"], ["c"]]
    assert len(split_sample_chunks(["a", "b"], 99)) == 2


def test_merge_scan_results_sums_switches():
    from flare_switch_qc import SwitchEvent, merge_scan_results, QcAccum

    e1 = SwitchEvent(
        sample="S1", hap=1, chrom="chr22", pos=100, ref=".", alt=".",
        anc_from=2, anc_to=3, gq=None, dp=None, rnc=".", prev_pos=50, bp_gap=50,
        is_flicker=True,
    )
    e2 = SwitchEvent(
        sample="S2", hap=1, chrom="chr22", pos=200, ref=".", alt=".",
        anc_from=3, anc_to=2, gq=None, dp=None, rnc=".", prev_pos=100, bp_gap=100,
    )
    meta1 = {
        "n_samples": 1, "sites_scanned": 10, "n_switches": 1,
        "n_flicker_switches": 1, "span_bp": 1000, "gq_threshold": 20, "dp_threshold": 5,
    }
    meta2 = {
        "n_samples": 1, "sites_scanned": 10, "n_switches": 1,
        "n_flicker_switches": 0, "span_bp": 1000, "gq_threshold": 20, "dp_threshold": 5,
    }
    events, qc, meta, tracts = merge_scan_results(
        [([e1], QcAccum(), meta1, []), ([e2], QcAccum(), meta2, [])]
    )
    assert len(events) == 2
    assert meta["n_samples"] == 2
    assert meta["n_switches"] == 2
    assert meta["n_flicker_switches"] == 1
    assert meta["jobs"] == 2
    assert qc.sites_scanned == 10


def test_format_id_header_lines_only():
    from flare_switch_qc import FORMAT_ID_RE

    assert FORMAT_ID_RE.match("##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"x\">").group(1) == "AN1"
    assert FORMAT_ID_RE.match("##contig=<ID=chr1,length=1>") is None
    # Substring in a description must not count as a FORMAT ID.
    assert FORMAT_ID_RE.match("##INFO=<ID=X,Description=\"##FORMAT=<ID=GQ, fake\">") is None


def test_sample_format_tags_requires_an():
    try:
        sample_format_tags({"GQ": True, "AN1": False, "AN2": True})
    except RuntimeError as exc:
        assert "AN1" in str(exc)
    else:
        raise AssertionError("expected RuntimeError")


def test_call_qc_filters_poor_switches():
    good = SwitchEvent(
        sample="S1",
        hap=1,
        chrom="chr22",
        pos=100,
        ref="A",
        alt="G",
        anc_from=2,
        anc_to=3,
        gq=40,
        dp=20,
        rnc="..",
        prev_pos=50,
        bp_gap=50,
    )
    bad_gq = SwitchEvent(
        sample="S1",
        hap=1,
        chrom="chr22",
        pos=200,
        ref="C",
        alt="T",
        anc_from=3,
        anc_to=2,
        gq=5,
        dp=20,
        rnc="..",
        prev_pos=100,
        bp_gap=100,
    )
    bad_rnc = SwitchEvent(
        sample="S2",
        hap=2,
        chrom="chr22",
        pos=300,
        ref="G",
        alt="A",
        anc_from=2,
        anc_to=1,
        gq=40,
        dp=20,
        rnc="II",
        prev_pos=200,
        bp_gap=100,
    )
    assert not switch_fails_call_qc(good, gq_threshold=20, dp_threshold=10)
    assert switch_fails_call_qc(bad_gq, gq_threshold=20, dp_threshold=10)
    assert switch_fails_call_qc(bad_rnc, gq_threshold=20, dp_threshold=10)
    passed, failed = filter_switches_call_qc(
        [good, bad_gq, bad_rnc], gq_threshold=20, dp_threshold=10
    )
    assert len(passed) == 1 and len(failed) == 2
    meta = {"n_samples": 2, "span_bp": 1_000_000, "gq_threshold": 20, "dp_threshold": 10}
    cqf = call_qc_filtered_tract_metrics(
        [good, bad_gq, bad_rnc], meta, gq_threshold=20, dp_threshold=10
    )
    assert cqf["n_switches_all"] == 3
    assert cqf["n_switches_fail_call_qc"] == 2
    assert cqf["frac_switches_fail_call_qc"] == 2 / 3
    assert abs(cqf["switches_per_hap_per_mb"] - 0.25) < 1e-9
    summary = summarize(QcAccum(), meta, [good, bad_gq, bad_rnc], [])
    assert summary["call_qc_filtered"]["n_switches_pass_call_qc"] == 1


if __name__ == "__main__":
    test_detects_switch_and_flicker()
    test_no_switch_when_anc_stable()
    test_summarize_enrichment()
    test_classify_gt_and_sentinel()
    test_summarize_het_enrichment_not_confounded_by_homref()
    test_nearest_het_for_homref_switch()
    test_nearest_het_is_current_when_switch_is_het()
    test_prev_het_does_not_cross_chromosomes()
    test_summarize_nearest_het_enrichment()
    test_homref_switch_records_incomplete_rnc()
    test_summarize_rnc_enrichment()
    test_expected_rate_from_t()
    test_build_tracts_includes_no_switch_haps_and_censoring()
    test_an_only_query_omits_gq_dp()
    test_split_sample_chunks()
    test_merge_scan_results_sums_switches()
    test_format_id_header_lines_only()
    test_sample_format_tags_requires_an()
    test_call_qc_filters_poor_switches()
    print("ok")
