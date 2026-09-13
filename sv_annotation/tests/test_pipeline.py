#!/usr/bin/env python3
"""End-to-end unit tests for the SV annotation script chain."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
FIX = Path(__file__).resolve().parent / "fixtures"


def run(*args: str) -> None:
    subprocess.check_call([sys.executable, *args])


@pytest.fixture
def work(tmp_path: Path) -> Path:
    return tmp_path


def test_ancestry_order():
    sys.path.insert(0, str(SCRIPTS))
    from sv_site_utils import order_samples_by_ancestry

    ordered = order_samples_by_ancestry({"S3": "afr", "S1": "eur", "S2": "amr"})
    assert ordered == ["S2", "S1", "S3"]


def test_binary_variant_magic(work: Path):
    sys.path.insert(0, str(SCRIPTS))
    from sv_site_utils import is_binary_variant_file

    gzip_path = work / "no_ext"
    gzip_path.write_bytes(b"\x1f\x8b\x08\x00")
    bcf_path = work / "also_no_ext"
    bcf_path.write_bytes(b"BCF\x02")
    vcf_path = FIX / "main.vcf"
    assert is_binary_variant_file(gzip_path)
    assert is_binary_variant_file(bcf_path)
    assert not is_binary_variant_file(vcf_path)


def test_pipeline_integrate_and_counts(work: Path):
    main_sites = work / "main.sites.tsv"
    bnd_sites = work / "bnd.sites.tsv"
    large_sites = work / "large.sites.tsv"
    carriers = work / "main.carriers.tsv"

    for src, vcf, out in (
        ("main", FIX / "main.vcf", main_sites),
        ("bnd", FIX / "bnd.vcf", bnd_sites),
        ("large", FIX / "large.vcf", large_sites),
    ):
        raw = work / f"{src}.raw.tsv"
        af = work / f"{src}.af.tsv"
        carr = work / f"{src}.carriers.tsv"
        run(
            str(SCRIPTS / "extract_sites_from_vcf.py"),
            "--vcf",
            str(vcf),
            "--source-vcf",
            src,
            "--phase",
            "phase2",
            "--out",
            str(raw),
        )
        run(
            str(SCRIPTS / "fill_af_and_carriers.py"),
            "--vcf",
            str(vcf),
            "--sites",
            str(raw),
            "--out-sites",
            str(af),
            "--out-carriers",
            str(carr),
        )
        run(
            str(SCRIPTS / "annotate_repeat_context.py"),
            "--sites",
            str(af),
            "--rmsk-bed",
            str(FIX / "rmsk.bed"),
            "--simple-repeat-bed",
            str(FIX / "simpleRepeat.bed"),
            "--segdup-bed",
            str(FIX / "segdup.bed"),
            "--cmrg-bed",
            str(FIX / "cmrg.bed"),
            "--out",
            str(work / f"{src}.region.tsv"),
        )
        if src == "main":
            run(
                str(SCRIPTS / "attach_caddsv_scores.py"),
                "--sites",
                str(work / f"{src}.region.tsv"),
                "--scores",
                str(FIX / "cadd_scores.tsv"),
                "--out",
                str(out),
            )
            # keep carriers from main
            carr.replace(carriers)
        else:
            (work / f"{src}.region.tsv").replace(out)

    unified = work / "unified.tsv"
    manifest = work / "manifest.json"
    run(
        str(SCRIPTS / "integrate_partitions.py"),
        "--main",
        str(main_sites),
        "--bnd",
        str(bnd_sites),
        "--large",
        str(large_sites),
        "--out",
        str(unified),
        "--manifest",
        str(manifest),
    )
    man = json.loads(manifest.read_text())
    assert man["n_main_suppressed_as_large_duplicate"] == 1
    assert man["n_output_sites"] == (
        man["n_main_input"] - 1 + man["n_bnd_input"] + man["n_large_input"]
    )

    text = unified.read_text()
    assert "overlap_main" not in text
    assert "overlap_large" in text
    assert "bnd1" in text

    counts = work / "counts.tsv"
    run(
        str(SCRIPTS / "manuscript_site_counts.py"),
        "--phase2-sites",
        str(unified),
        "--out-tsv",
        str(counts),
        "--out-json",
        str(work / "counts.json"),
    )
    lines = {
        r.split("\t")[0]: r.split("\t")
        for r in counts.read_text().splitlines()[1:]
    }
    assert lines["breakends"][2] == "1"
    assert lines["large_events_gt10kb"][2] == "2"
    assert lines["deletions"][2].startswith("1;")  # only del1 (overlap suppressed)
    assert lines["deletions_us"][2] == "0; 0"
    assert lines["deletions_rm"][2].startswith("1;")
    assert lines["deletions_sd"][2] == "0; 0"
    assert lines["deletions_sr"][2] == "0; 0"
    assert lines["deletions_cmrg"][2].startswith("1;")
    assert lines["insertions"][2].startswith("1;")
    assert lines["insertions_us"][2] == "0; 0"
    assert lines["insertions_rm"][2] == "0; 0"
    assert lines["insertions_sd"][2] == "0; 0"
    assert lines["insertions_sr"][2].startswith("1;")
    assert lines["insertions_cmrg"][2] == "0; 0"
    assert lines["deletions_context_pct"][2] == "0.0 / 100.0 / 0.0 / 0.0 / 100.0"
    assert lines["insertions_context_pct"][2] == "0.0 / 0.0 / 0.0 / 100.0 / 0.0"
    assert lines["duplications"][2].startswith("1;")
    assert lines["inversions"][2].startswith("1;")

    disc = work / "discovery.tsv"
    run(
        str(SCRIPTS / "ebert_discovery.py"),
        "--sites",
        str(main_sites),
        "--carriers",
        str(carriers),
        "--sample-ancestry",
        str(FIX / "sample_ancestry.tsv"),
        "--strata",
        "none",
        "--include-sources",
        "main",
        "--out",
        str(disc),
    )
    # First sample in ancestry order should be amr (S2)
    first = disc.read_text().splitlines()[1].split("\t")
    assert first[1] == "S2"
    assert first[2] == "amr"

    # Multi-strata one-pass matches single-stratum output
    run(
        str(SCRIPTS / "ebert_discovery.py"),
        "--sites",
        str(main_sites),
        "--carriers",
        str(carriers),
        "--sample-ancestry",
        str(FIX / "sample_ancestry.tsv"),
        "--strata-list",
        "none,region,cadd",
        "--include-sources",
        "main",
        "--out-prefix",
        str(work / "multi"),
    )
    assert (work / "multi.discovery.tsv").read_text() == disc.read_text()
    assert (work / "multi.discovery.region.tsv").exists()
    assert (work / "multi.discovery.cadd.tsv").exists()


def test_cadd_bins():
    sys.path.insert(0, str(SCRIPTS))
    from sv_site_utils import cadd_sv_bin

    assert cadd_sv_bin(None) == "unscored"
    assert cadd_sv_bin(9.9) == "low"
    assert cadd_sv_bin(10) == "mid"
    assert cadd_sv_bin(20) == "high"


def test_caddsv_drops_megabase_dup_inv(work: Path):
    sys.path.insert(0, str(SCRIPTS))
    from run_caddsv import filter_caddsv_bed

    src = work / "in.bed"
    dst = work / "out.bed"
    src.write_text(
        "chr1\t0\t1000\tDEL\n"
        "chr1\t0\t2000000\tDUP\n"
        "chr1\t0\t5000\tINV\n"
        "chr1\t0\t5000000\tINV\n"
        "chr2\t10\t80\tINS\n"
    )
    kept = filter_caddsv_bed(src, dst, max_dup_inv_bp=1_000_000)
    rows = [line.split("\t")[3] for line in dst.read_text().splitlines()]
    assert kept == 3
    assert rows == ["DEL", "INV", "INS"]


def test_caddsv_input_is_v2_coordinate_bed(work: Path):
    raw = work / "raw.tsv"
    run(
        str(SCRIPTS / "extract_sites_from_vcf.py"),
        "--vcf",
        str(FIX / "main.vcf"),
        "--source-vcf",
        "main",
        "--phase",
        "phase2",
        "--out",
        str(raw),
    )
    bed = work / "caddsv.bed"
    run(
        str(SCRIPTS / "prepare_caddsv_input.py"),
        "--sites",
        str(raw),
        "--out-bed",
        str(bed),
    )
    rows = [line.split("\t") for line in bed.read_text().splitlines()]
    assert rows
    assert all(len(row) == 4 for row in rows)
    assert {row[3] for row in rows} == {"DEL", "DUP", "INS", "INV"}
    assert "ins1" not in bed.read_text()


def test_integrate_main_only_without_bnd_or_large(work: Path):
    """Phase 1 has no bnd/large partitions; integrate must not require them."""
    raw = work / "main.raw.tsv"
    main_sites = work / "main.sites.tsv"
    carr = work / "main.carriers.tsv"
    run(
        str(SCRIPTS / "extract_sites_from_vcf.py"),
        "--vcf",
        str(FIX / "main.vcf"),
        "--source-vcf",
        "main",
        "--phase",
        "phase1",
        "--out",
        str(raw),
    )
    run(
        str(SCRIPTS / "fill_af_and_carriers.py"),
        "--vcf",
        str(FIX / "main.vcf"),
        "--sites",
        str(raw),
        "--out-sites",
        str(main_sites),
        "--out-carriers",
        str(carr),
    )
    unified = work / "unified.tsv"
    manifest = work / "manifest.json"
    run(
        str(SCRIPTS / "integrate_partitions.py"),
        "--main",
        str(main_sites),
        "--out",
        str(unified),
        "--manifest",
        str(manifest),
    )
    man = json.loads(manifest.read_text())
    assert man["n_bnd_input"] == 0
    assert man["n_large_input"] == 0
    assert man["n_main_suppressed_as_large_duplicate"] == 0
    assert man["n_output_sites"] == man["n_main_input"]
    assert man["n_output_sites"] > 0


def test_phase1_counts_mark_bnd_unavailable(work: Path):
    raw = work / "p1.raw.tsv"
    run(
        str(SCRIPTS / "extract_sites_from_vcf.py"),
        "--vcf",
        str(FIX / "main.vcf"),
        "--source-vcf",
        "main",
        "--phase",
        "phase1",
        "--out",
        str(raw),
    )
    # minimal af pass
    af = work / "p1.af.tsv"
    carr = work / "p1.carr.tsv"
    run(
        str(SCRIPTS / "fill_af_and_carriers.py"),
        "--vcf",
        str(FIX / "main.vcf"),
        "--sites",
        str(raw),
        "--out-sites",
        str(af),
        "--out-carriers",
        str(carr),
    )
    counts = work / "p1.counts.tsv"
    run(
        str(SCRIPTS / "manuscript_site_counts.py"),
        "--phase1-sites",
        str(af),
        "--out-tsv",
        str(counts),
        "--out-json",
        str(work / "p1.json"),
    )
    rows = {
        r.split("\t")[0]: r.split("\t")
        for r in counts.read_text().splitlines()[1:]
    }
    assert rows["breakends"][1] == "—"
    assert rows["large_events_gt10kb"][1] == "—"
    assert rows["deletions"][1] != "—"
