#!/usr/bin/env python3
"""Unit tests for flare_switch_qc (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_switch_qc import (  # noqa: E402
    SampleTracker,
    SwitchEvent,
    mark_flickers,
    maybe_record_switch,
    summarize,
    QcAccum,
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


if __name__ == "__main__":
    test_detects_switch_and_flicker()
    test_no_switch_when_anc_stable()
    test_summarize_enrichment()
    print("ok")
