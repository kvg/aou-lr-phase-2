#!/usr/bin/env python3
"""Tests for trgt_defs_from_vcf.py.

The fixture is a TRGT-shaped VCF built to exercise the real bcftools path,
including a compound locus whose TRID and MOTIFS are comma-joined and whose
STRUC is the symbolic form TRExplorer uses.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import trgt_defs_from_vcf as tdv  # noqa: E402

SCRIPT = Path(__file__).resolve().parent / "trgt_defs_from_vcf.py"

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=TRID,Number=1,Type=String,Description="Locus id">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">
##INFO=<ID=STRUC,Number=1,Type=String,Description="Structure">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""

# (POS, INFO, ) — mirrors the real file: a simple locus, a compound locus, and
# one record with STRUC missing.
ROWS = [
    (66431, "TRID=1-66431-66446-TATAT;END=66446;MOTIFS=TATAT;STRUC=<TR22>"),
    (66160, "TRID=1-66160-66175-TA,1-66192-66201-TA;END=66415;MOTIFS=ATATA,TATAT,TA;STRUC=<VC2>"),
    (70351, "TRID=1-70351-70360-A;END=70360;MOTIFS=A;STRUC=."),
]


@pytest.fixture
def vcf(tmp_path: Path) -> Path:
    p = tmp_path / "trgt.vcf"
    lines = [HEADER.rstrip("\n")]
    for pos, info in ROWS:
        lines.append("\t".join(("chr1", str(pos), ".", "A" * 10, "<TR>", ".", ".", info,
                                "GT", "0/1")))
    p.write_text("\n".join(lines) + "\n")
    return p


def run_script(vcf: Path, out: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--vcf", str(vcf), "--out", str(out), *extra],
        capture_output=True, text=True)


def test_definition_round_trips_pos_to_bed_start():
    out = tdv.definition("chr20", "10000930", "10000950", "20-10000930-10000950-A", "A", "<TR1>")
    fields = out.split("\t")
    assert fields[:3] == ["chr20", "10000929", "10000950"], "BED start is POS-1, end is INFO END"
    assert fields[3] == "ID=20-10000930-10000950-A;MOTIFS=A;STRUC=<TR1>"
    assert len(fields) == 4, "trgt --repeats accepts exactly four fields"


def test_output_passes_the_workflow_guard(vcf: Path, tmp_path: Path):
    """The WDL refuses a first row that is not 4 fields containing ID=."""
    out = tmp_path / "defs.bed"
    proc = run_script(vcf, out)
    assert proc.returncode == 0, proc.stderr
    first = out.read_text().splitlines()[0]
    assert len(first.split("\t")) == 4 and "ID=" in first


def test_compound_locus_keeps_its_comma_joined_trid_and_motifs(vcf: Path, tmp_path: Path):
    out = tmp_path / "defs.bed"
    assert run_script(vcf, out).returncode == 0
    rows = {l.split("\t")[1]: l.split("\t")[3] for l in out.read_text().splitlines()}
    compound = rows["66159"]
    assert compound == ("ID=1-66160-66175-TA,1-66192-66201-TA;"
                        "MOTIFS=ATATA,TATAT,TA;STRUC=<VC2>")


def test_record_missing_struc_is_skipped_and_counted(vcf: Path, tmp_path: Path):
    out = tmp_path / "defs.bed"
    proc = run_script(vcf, out)
    assert proc.returncode == 0
    assert len(out.read_text().splitlines()) == 2, "the STRUC-less record is dropped"
    assert "skipped 1 record(s)" in proc.stderr
    assert "70351" in proc.stderr


def test_gz_output_is_compressed(vcf: Path, tmp_path: Path):
    out = tmp_path / "defs.bed.gz"
    assert run_script(vcf, out).returncode == 0
    with gzip.open(out, "rt") as fh:
        assert fh.readline().startswith("chr1\t")


def test_expect_loci_mismatch_is_fatal(vcf: Path, tmp_path: Path):
    proc = run_script(vcf, tmp_path / "defs.bed", "--expect-loci", "14618")
    assert proc.returncode != 0
    assert "--expect-loci" in proc.stderr
    assert "only the loci TRGT emitted" in proc.stderr


def test_expect_loci_match_passes(vcf: Path, tmp_path: Path):
    proc = run_script(vcf, tmp_path / "defs.bed", "--expect-loci", "2")
    assert proc.returncode == 0, proc.stderr


def test_a_vcf_with_no_struc_at_all_is_fatal(tmp_path: Path):
    p = tmp_path / "nostruc.vcf"
    p.write_text(HEADER + "\t".join((
        "chr1", "100", ".", "AAAA", "<TR>", ".", ".",
        "TRID=x;END=110;MOTIFS=A;STRUC=.", "GT", "0/1")) + "\n")
    proc = run_script(p, tmp_path / "defs.bed")
    assert proc.returncode != 0
    assert "not a TRGT output VCF" in proc.stderr


def test_bcftools_failure_is_surfaced(tmp_path: Path):
    with pytest.raises(SystemExit) as exc:
        list(tdv.query_lines(str(tmp_path / "missing.vcf.gz")))
    assert "bcftools query failed" in str(exc.value)
