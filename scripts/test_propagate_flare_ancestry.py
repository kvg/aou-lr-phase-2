#!/usr/bin/env python3
"""Unit tests for propagate_flare_ancestry (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from propagate_flare_ancestry import (  # noqa: E402
    Region,
    flare_query_region,
    main,
    propagate_ancestry,
)


def _records(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln and not ln.startswith("#")]


def _write_vcf(path: Path, samples: str, records: str, extra_header: str = "") -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr22>\n"
        + extra_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + samples
        + "\n"
        + records
    )


def test_felixla_intervals_and_preserves_target_gt(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1\tS2",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\t0|0:2:2\n"
        "chr22\t200\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t1|1:3:3\t0|1:2:3\n",
        extra_header=(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="A">\n'
            '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="B">\n'
        ),
    )
    _write_vcf(
        target,
        "S1\tS2",
        "chr22\t50\tsv_before\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0|1\t0|0\n"
        "chr22\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n"
        "chr22\t150\tinterstitial\tG\tA\t.\tPASS\t.\tGT\t1|0\t0|1\n"
        "chr22\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|1\n"
        "chr22\t250\tafter\tT\tC\t.\tPASS\t.\tGT\t0|0\t1|0\n",
    )
    out = tmp_path / "out.vcf"
    stats = propagate_ancestry(flare, target, out)
    assert stats.target_sites == 5
    assert stats.annotated_sites == 5
    assert stats.sites_missing_ancestry == 0
    assert stats.exact_pos_overlap == 2
    assert stats.sites_before_first_flare == 1
    assert stats.sites_after_last_flare == 1
    recs = _records(out.read_text())
    assert len(recs) == 5
    fmt = recs[0].split("\t")[8]
    assert fmt == "GT:AN1:AN2"
    # FELIXla: first marker covers 1..100
    assert recs[0].split("\t")[9] == "0|1:2:3"
    assert recs[0].split("\t")[7].endswith("FLARE_POS=100") or "FLARE_POS=100" in recs[0]
    # exact match at 100
    assert recs[1].split("\t")[9] == "0|1:2:3"
    # interstitial 150 uses marker at 200, not left-carry from 100
    assert recs[2].split("\t")[9] == "1|0:3:3"
    assert "FLARE_POS=200" in recs[2]
    assert recs[3].split("\t")[9] == "1|1:3:3"
    # last state extends past 200
    assert recs[4].split("\t")[9] == "0|0:3:3"
    assert recs[4].split("\t")[10] == "1|0:2:3"


def test_sample_name_alignment_not_column_order(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1\tS2",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\t0|0:4:5\n",
    )
    _write_vcf(
        target,
        "S2\tS1",
        "chr22\t150\t.\tG\tA\t.\t.\t.\tGT\t1|0\t0|1\n",
    )
    out = tmp_path / "out.vcf"
    stats = propagate_ancestry(flare, target, out)
    assert stats.n_shared_samples == 2
    recs = _records(out.read_text())
    s2, s1 = recs[0].split("\t")[9:11]
    assert s2 == "1|0:4:5"
    assert s1 == "0|1:2:3"


def test_missing_flare_sample_and_missing_contig(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n",
    )
    _write_vcf(
        target,
        "S1\tS_extra",
        "chr22\t150\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|0\n",
    )
    out = tmp_path / "out.vcf"
    stats = propagate_ancestry(flare, target, out)
    assert stats.n_shared_samples == 1
    recs = _records(out.read_text())
    assert recs[0].split("\t")[9] == "0|1:2:3"
    assert recs[0].split("\t")[10] == "1|0:.:."

    empty_flare = tmp_path / "flare_chr1.vcf"
    _write_vcf(
        empty_flare,
        "S1",
        "chr1\t10\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|0:1:1\n",
    )
    missing = tmp_path / "missing.tsv"
    out2 = tmp_path / "out2.vcf"
    stats2 = propagate_ancestry(empty_flare, target, out2, missing_tsv=missing)
    assert stats2.annotated_sites == 0
    assert stats2.sites_missing_ancestry == 1
    assert missing.read_text().splitlines()[1] == "chr22\t150\tG\tA"
    assert _records(out2.read_text())[0].split("\t")[9] == "0|1:.:."


def test_region_keeps_only_interval(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n"
        "chr22\t200\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:4:5\n"
        "chr22\t300\t.\tT\tC\t.\t.\t.\tGT:AN1:AN2\t0|0:6:7\n",
    )
    _write_vcf(
        target,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\n"
        "chr22\t150\t.\tG\tA\t.\t.\t.\tGT\t1|0\n"
        "chr22\t250\t.\tC\tT\t.\t.\t.\tGT\t0|0\n",
    )
    out = tmp_path / "out.vcf"
    stats = propagate_ancestry(
        flare, target, out, region=Region.parse("chr22:140-200")
    )
    recs = _records(out.read_text())
    assert [int(r.split("\t")[1]) for r in recs] == [150]
    assert recs[0].split("\t")[9] == "1|0:4:5"
    assert stats.target_sites == 1


def test_min_coverage_fails_after_writing(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1",
        "chr1\t10\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|0:1:1\n",
    )
    _write_vcf(
        target,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\n",
    )
    out = tmp_path / "out.vcf"
    stats_json = tmp_path / "stats.json"
    with pytest.raises(SystemExit, match="coverage too low"):
        main(
            [
                "--flare",
                str(flare),
                "--target",
                str(target),
                "--output",
                str(out),
                "--stats-json",
                str(stats_json),
                "--min-coverage-rate",
                "0.95",
            ]
        )
    assert out.is_file()
    assert '"annotated_sites": 0' in stats_json.read_text()


def test_unsorted_flare_fails(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1",
        "chr22\t200\t.\tC\tT\t.\t.\t.\tGT:AN1:AN2\t0|0:2:2\n"
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n",
    )
    _write_vcf(
        target,
        "S1",
        "chr22\t150\t.\tG\tA\t.\t.\t.\tGT\t0|1\n",
    )
    with pytest.raises(ValueError, match="unsorted"):
        propagate_ancestry(flare, target, tmp_path / "out.vcf")


def test_cli_writes_stats(tmp_path: Path):
    flare = tmp_path / "flare.vcf"
    target = tmp_path / "target.vcf"
    _write_vcf(
        flare,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT:AN1:AN2\t0|1:2:3\n",
    )
    _write_vcf(
        target,
        "S1",
        "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1\n"
        "chr22\t150\t.\tG\tA\t.\t.\t.\tGT\t0|0\n",
    )
    out = tmp_path / "out.vcf"
    stats_json = tmp_path / "stats.json"
    rc = main(
        [
            "--flare",
            str(flare),
            "--target",
            str(target),
            "--output",
            str(out),
            "--stats-json",
            str(stats_json),
        ]
    )
    assert rc == 0
    text = stats_json.read_text()
    assert '"annotated_sites": 2' in text
    assert '"exact_pos_overlap": 1' in text


def test_flare_query_region_left_pads():
    assert flare_query_region("chr22") == "chr22"
    assert flare_query_region("chr22:16000000-17000000") == "chr22:1-17000000"
    assert Region.parse("chr22:16000000-17000000").contains("chr22", 16000000)
    assert not Region.parse("chr22:16000000-17000000").contains("chr22", 15999999)
