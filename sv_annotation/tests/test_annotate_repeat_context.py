"""Breakpoint context labels (Zhao / Xuefang), not whole-interval unions."""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from annotate_repeat_context import (  # noqa: E402
    contains_point,
    context_label,
    covered_bases,
    load_intervals,
    reference_interval,
)
from sv_site_utils import site_row, write_site_header  # noqa: E402


def _iv(*spans: tuple[int, int]) -> list[tuple[int, int]]:
    return list(spans)


def test_breakpoint_sr_ignores_rm_in_the_body():
    tracks = {
        "RM": _iv((102, 108)),
        "SD": [],
        "SR": _iv((99, 100)),
    }
    label, hits = context_label("DEL", 99, 110, tracks)
    assert label == "SR"
    assert hits == {"RM": False, "SD": False, "SR": True}


def test_sr_overrides_rm_when_both_breakpoints_hit():
    tracks = {
        "RM": _iv((90, 120)),
        "SD": [],
        "SR": _iv((90, 100)),
    }
    label, hits = context_label("INS", 99, 100, tracks)
    assert label == "SR"
    assert hits["RM"] and hits["SR"]


def test_long_deletion_uses_body_coverage_not_breakpoints():
    # Span (999, 7000) is 6001 bp. SD covers 3501 bp (>0.5). RM hits only the
    # left breakpoint, which the body rule ignores.
    tracks = {
        "RM": _iv((999, 1000)),
        "SD": _iv((999, 4500)),
        "SR": [],
    }
    label, hits = context_label("DEL", 999, 7000, tracks)
    assert label == "SD"
    assert hits == {"RM": False, "SD": True, "SR": False}


def test_span_of_5000_stays_on_breakpoints():
    # Exactly 5000 bp is not the body rule (threshold is span > 5000).
    # RM covers the middle; SR hits the left breakpoint.
    tracks = {
        "RM": _iv((100, 4000)),
        "SD": [],
        "SR": _iv((0, 1)),
    }
    label, _hits = context_label("DEL", 0, 5000, tracks)
    assert label == "SR"


def test_long_insertion_stays_on_the_anchor():
    site = {
        "pos": 1000,
        "end": 8000,
        "svtype": "INS",
        "svlen": 7000,
    }
    start0, end0 = reference_interval(site)
    assert (start0, end0) == (999, 1000)
    tracks = {
        "RM": [],
        "SD": _iv((2000, 8000)),
        "SR": _iv((999, 1000)),
    }
    label, hits = context_label("INS", start0, end0, tracks)
    assert label == "SR"
    assert hits["SD"] is False


def test_overlapping_copies_do_not_inflate_body_fraction(tmp_path: Path):
    bed = tmp_path / "rm.bed"
    # Two identical copies of a 4000 bp interval over a 10000 bp deletion.
    bed.write_text("chr1\t0\t4000\nchr1\t0\t4000\n")
    merged = load_intervals(str(bed))["chr1"]
    assert covered_bases(merged, 0, 10000) == 4000
    label, hits = context_label("DEL", 0, 10000, {"RM": merged, "SD": [], "SR": []})
    assert label == "US"
    assert hits["RM"] is False


def test_point_on_interval_end_is_outside():
    assert contains_point([(10, 20)], 10)
    assert not contains_point([(10, 20)], 20)


def test_cli_writes_exclusive_class_and_cmrg(tmp_path: Path):
    sites = tmp_path / "sites.tsv"
    with sites.open("w") as fh:
        write_site_header(fh)
        fh.write(
            site_row(
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "end": 200,
                    "id": "del",
                    "svtype": "DEL",
                    "svlen": -100,
                    "region_class": "US",
                }
            )
        )
    rmsk = tmp_path / "rmsk.bed"
    sr = tmp_path / "sr.bed"
    sd = tmp_path / "sd.bed"
    cmrg = tmp_path / "cmrg.bed"
    rmsk.write_text("chr1\t99\t100\n")
    sr.write_text("chr1\t99\t100\n")
    sd.write_text("")
    cmrg.write_text("chr1\t150\t180\n")
    out = tmp_path / "out.tsv"
    from annotate_repeat_context import main
    import sys

    argv = sys.argv
    sys.argv = [
        "annotate_repeat_context.py",
        "--sites",
        str(sites),
        "--rmsk-bed",
        str(rmsk),
        "--simple-repeat-bed",
        str(sr),
        "--segdup-bed",
        str(sd),
        "--cmrg-bed",
        str(cmrg),
        "--out",
        str(out),
    ]
    try:
        main()
    finally:
        sys.argv = argv
    row = out.read_text().splitlines()[1].split("\t")
    header = out.read_text().splitlines()[0].split("\t")
    got = dict(zip(header, row))
    assert got["region_class"] == "SR"
    assert got["hit_rmsk"] == "true"
    assert got["hit_simpleRepeat"] == "true"
    assert got["hit_genomicSuperDups"] == "false"
    assert got["hit_cmrg"] == "true"
