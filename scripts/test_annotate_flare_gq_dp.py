#!/usr/bin/env python3
"""Unit tests for annotate_flare_gq_dp (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
TESTDATA = REPO / "propagate_annotations" / "testdata"
sys.path.insert(0, str(SCRIPTS))

from annotate_flare_gq_dp import (  # noqa: E402
    Region,
    annotate_flare,
    extract_tags,
    main,
    match_dv_record,
)


def _records(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln and not ln.startswith("#")]


def test_extract_tags_from_glnexus_format():
    fmt = "GT:DP:AD:GQ:PL:RNC"
    sample = "0/1:20:10,10:30:0,30,60:.."
    assert extract_tags(fmt, sample, ("GQ", "DP")) == {"GQ": "30", "DP": "20"}
    assert extract_tags(fmt, sample, ("GQ", "DP", "RNC")) == {
        "GQ": "30",
        "DP": "20",
        "RNC": "..",
    }
    missing = "0/0:0:0,0:1:0,0,0:II"
    assert extract_tags(fmt, missing, ("RNC",)) == {"RNC": "II"}


def test_copies_gq_dp_and_keeps_flare_gt_an(tmp_path: Path):
    out = tmp_path / "out.vcf"
    stats = annotate_flare(
        TESTDATA / "flare.vcf",
        TESTDATA / "snv_indel.vcf",
        out,
    )
    assert stats.flare_sites == 3
    assert stats.matched_sites == 3
    assert stats.unmatched_flare_sites == 0
    assert stats.match_rate == 1.0
    assert stats.n_shared_samples == 2

    text = out.read_text()
    assert "##FORMAT=<ID=GQ," in text
    assert "##FORMAT=<ID=RNC," in text
    assert "##FORMAT=<ID=AN1," in text
    recs = _records(text)
    assert len(recs) == 3
    fmt = recs[0].split("\t")[8]
    assert fmt == "GT:AN1:AN2:GQ:DP:RNC"
    s1 = recs[0].split("\t")[9]
    # FLARE GT/AN preserved; GQ=30 DP=20 RNC=.. from DeepVariant S1 at pos 100
    assert s1 == "0|1:2:3:30:20:.."
    s2 = recs[0].split("\t")[10]
    assert s2 == "0|0:2:2:25:15:.."
    # pos 300 S1 is GLnexus incomplete (RNC=II)
    assert recs[2].split("\t")[9] == "0|0:2:2:5:2:II"


def test_unmatched_flare_site_gets_missing_gq_dp(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    dv = tmp_path / "dv.vcf"
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"A\">\n"
        "##FORMAT=<ID=AN2,Number=1,Type=Integer,Description=\"B\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n"
        "chr22\t150\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:2:2\n"
    )
    dv.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\trs1\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:20:30\n"
        "chr22\t200\trs2\tG\tA\t.\t.\t.\tGT:DP:GQ\t0/0:10:15\n"
    )
    out = tmp_path / "out.vcf"
    stats = annotate_flare(flare, dv, out)
    assert stats.matched_sites == 1
    assert stats.unmatched_flare_sites == 1
    recs = _records(out.read_text())
    assert recs[1].split("\t")[9] == "0|0:2:2:.:.:."


def test_sample_name_alignment_not_column_order(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    dv = tmp_path / "dv.vcf"
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\t0|0:2:2\n"
    )
    dv.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/0:15:25\t0/1:20:30\n"
    )
    out = tmp_path / "out.vcf"
    stats = annotate_flare(flare, dv, out)
    assert stats.n_shared_samples == 2
    recs = _records(out.read_text())
    s1, s2 = recs[0].split("\t")[9:11]
    assert s1.endswith(":30:20:.")
    assert s2.endswith(":25:15:.")


def test_region_keeps_only_interval(tmp_path: Path):
    out = tmp_path / "out.vcf"
    stats = annotate_flare(
        TESTDATA / "flare.vcf",
        TESTDATA / "snv_indel.vcf",
        out,
        region=Region.parse("chr22:150-250"),
    )
    recs = _records(out.read_text())
    assert [int(r.split("\t")[1]) for r in recs] == [200]
    assert stats.flare_sites == 1
    assert stats.matched_sites == 1


def test_min_match_rate_fails_after_writing(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    dv = tmp_path / "dv.vcf"
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n"
        "chr22\t200\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:2:2\n"
    )
    dv.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:20:30\n"
    )
    out = tmp_path / "out.vcf"
    stats_json = tmp_path / "stats.json"
    with pytest.raises(SystemExit, match="match rate too low"):
        main(
            [
                "--flare",
                str(flare),
                "--deepvariant",
                str(dv),
                "--output",
                str(out),
                "--stats-json",
                str(stats_json),
                "--min-match-rate",
                "0.95",
            ]
        )
    assert out.is_file()
    assert '"matched_sites": 1' in stats_json.read_text()


def test_unsorted_deepvariant_fails(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    dv = tmp_path / "dv.vcf"
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n"
        "chr22\t300\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:2:2\n"
    )
    dv.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t300\t.\tC\tT\t.\t.\t.\tGT:DP:GQ\t0/0:10:15\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:20:30\n"
    )
    with pytest.raises(ValueError, match="unsorted"):
        annotate_flare(flare, dv, tmp_path / "out.vcf")


def test_cli_writes_stats(tmp_path: Path):
    out = tmp_path / "out.vcf"
    stats_json = tmp_path / "stats.json"
    rc = main(
        [
            "--flare",
            str(TESTDATA / "flare.vcf"),
            "--deepvariant",
            str(TESTDATA / "snv_indel.vcf"),
            "--output",
            str(out),
            "--stats-json",
            str(stats_json),
        ]
    )
    assert rc == 0
    assert out.is_file()
    assert '"matched_sites": 3' in stats_json.read_text()


def test_cli_output_dash_writes_stdout(tmp_path: Path, capsys):
    stats_json = tmp_path / "stats.json"
    rc = main(
        [
            "--flare",
            str(TESTDATA / "flare.vcf"),
            "--deepvariant",
            str(TESTDATA / "snv_indel.vcf"),
            "--output",
            "-",
            "--stats-json",
            str(stats_json),
        ]
    )
    assert rc == 0
    text = capsys.readouterr().out
    assert text.startswith("##fileformat=VCFv4.2")
    assert "##FORMAT=<ID=GQ," in text
    assert len(_records(text)) == 3
    assert '"matched_sites": 3' in stats_json.read_text()


def test_region_parse():
    r = Region.parse("chr22:16000000-17000000")
    assert r.contains("chr22", 16000000)
    assert r.contains("chr22", 17000000)
    assert not r.contains("chr22", 15999999)
    assert Region.parse("chr22").contains("chr22", 1)


def test_multiallelic_dv_alt_matches_split_flare(tmp_path: Path):
    rec, how = match_dv_record(
        "A",
        "G",
        [["chr22", "100", ".", "A", "G,T", ".", ".", ".", "GT:GQ:DP", "0/1:30:20"]],
    )
    assert how == "multiallelic"
    assert rec is not None
    flare = tmp_path / "flare.vcf"
    dv = tmp_path / "dv.vcf"
    unmatched = tmp_path / "unmatched.tsv"
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n"
        "chr22\t200\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:2:2\n"
    )
    dv.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr22\t100\t.\tA\tG,T\t.\t.\t.\tGT:DP:GQ\t0/1:20:30\n"
    )
    out = tmp_path / "out.vcf"
    stats = annotate_flare(flare, dv, out, unmatched_tsv=unmatched)
    assert stats.matched_sites == 1
    assert stats.matched_multiallelic == 1
    assert stats.unmatched_flare_sites == 1
    assert unmatched.read_text().splitlines()[1] == "chr22\t200\tC\tT"
    assert _records(out.read_text())[0].split("\t")[9] == "0|1:2:3:30:20:."
