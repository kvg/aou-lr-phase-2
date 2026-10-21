#!/usr/bin/env python3
"""Tests for make_trgt_phase_table.py."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_trgt_phase_table as mpt  # noqa: E402

SCRIPT = Path(__file__).resolve().parent / "make_trgt_phase_table.py"

BAMS_HEADER = "entity:aou2_v1_phased_bams_id\tgcs_bai_path\tgcs_bam_path\tsample\n"


def write_bams(tmp_path: Path, rows: list[tuple[str, str, str]]) -> Path:
    p = tmp_path / "bams.tsv"
    body = "".join(f"{i}\t{bai}\t{bam}\t{i}\n" for i, bai, bam in rows)
    p.write_text(BAMS_HEADER + body)
    return p


def write_trgt(tmp_path: Path, rows: list[tuple[str, str]]) -> Path:
    p = tmp_path / "trgt.txt"
    p.write_text("".join(f"{i}\t{uri}\n" for i, uri in rows))
    return p


def run_script(tmp_path: Path, *extra: str) -> tuple[str, list[str]]:
    out = tmp_path / "trgt_phase.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--bams-table", str(tmp_path / "bams.tsv"),
         "--trgt-table", str(tmp_path / "trgt.txt"), "--out", str(out), *extra],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout + proc.stderr, out.read_text().splitlines()


def test_only_the_intersection_is_written(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", "gs://b/a.bam"),
                          ("b", "gs://b/b.bam.bai", "gs://b/b.bam")])
    write_trgt(tmp_path, [("b", "gs://t/b.vcf.gz"), ("c", "gs://t/c.vcf.gz")])
    stdout, lines = run_script(tmp_path)
    assert [l.split("\t")[0] for l in lines[1:]] == ["b"]
    assert "ROWS WRITTEN (both sources):    1" in stdout
    assert "BAM but no TRGT call:         1" in stdout
    assert "TRGT call but no usable BAM:  1" in stdout


def test_header_and_row_shape(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", "gs://b/a.bam")])
    write_trgt(tmp_path, [("a", "gs://t/a.vcf.gz")])
    _, lines = run_script(tmp_path, "--entity-type", "trgt_phase")
    assert lines[0] == ("entity:trgt_phase_id\tsample\tbam\tbam_index\ttrgt_vcf\ttrgt_vcf_index")
    assert lines[1].split("\t") == [
        "a", "a", "gs://b/a.bam", "gs://b/a.bam.bai", "gs://t/a.vcf.gz", "gs://t/a.vcf.gz.tbi"]


def test_rows_follow_the_bam_table_order(tmp_path: Path):
    write_bams(tmp_path, [("c", "gs://b/c.bam.bai", "gs://b/c.bam"),
                          ("a", "gs://b/a.bam.bai", "gs://b/a.bam"),
                          ("b", "gs://b/b.bam.bai", "gs://b/b.bam")])
    write_trgt(tmp_path, [("a", "gs://t/a.vcf.gz"), ("b", "gs://t/b.vcf.gz"),
                          ("c", "gs://t/c.vcf.gz")])
    _, lines = run_script(tmp_path)
    assert [l.split("\t")[0] for l in lines[1:]] == ["c", "a", "b"]


def test_blank_bam_cell_excludes_the_sample(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", ""),
                          ("b", "", "gs://b/b.bam"),
                          ("c", "gs://b/c.bam.bai", "gs://b/c.bam")])
    write_trgt(tmp_path, [("a", "gs://t/a.vcf.gz"), ("b", "gs://t/b.vcf.gz"),
                          ("c", "gs://t/c.vcf.gz")])
    stdout, lines = run_script(tmp_path)
    assert [l.split("\t")[0] for l in lines[1:]] == ["c"]
    assert "missing a BAM or index cell:  2" in stdout


def test_emitted_duplicate_fails_until_resolved(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", "gs://b/a.bam")])
    write_trgt(tmp_path, [("a", "gs://t/early.vcf.gz"), ("a", "gs://t/late.vcf.gz")])
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--bams-table", str(tmp_path / "bams.tsv"),
         "--trgt-table", str(tmp_path / "trgt.txt"), "--out", str(tmp_path / "o.tsv")],
        capture_output=True, text=True)
    assert proc.returncode != 0
    assert "gs://t/early.vcf.gz" in proc.stderr and "--on-duplicate" in proc.stderr

    stdout, lines = run_script(tmp_path, "--on-duplicate", "last")
    assert lines[1].split("\t")[4] == "gs://t/late.vcf.gz"


def test_duplicate_outside_the_intersection_does_not_block(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", "gs://b/a.bam")])
    write_trgt(tmp_path, [("a", "gs://t/a.vcf.gz"),
                          ("z", "gs://t/early.vcf.gz"), ("z", "gs://t/late.vcf.gz")])
    stdout, lines = run_script(tmp_path)
    assert [l.split("\t")[0] for l in lines[1:]] == ["a"]
    assert "ids with >1 TRGT path:        1" in stdout


def test_missing_bam_column_is_fatal_and_lists_columns(tmp_path: Path):
    write_bams(tmp_path, [("a", "gs://b/a.bam.bai", "gs://b/a.bam")])
    with pytest.raises(SystemExit) as exc:
        mpt.read_bams_table(str(tmp_path / "bams.tsv"), "nope", "gcs_bai_path")
    assert "no column 'nope'" in str(exc.value)
    assert "gcs_bam_path" in str(exc.value)


def test_non_terra_header_is_rejected(tmp_path: Path):
    p = tmp_path / "bams.tsv"
    p.write_text("sample\tgcs_bam_path\tgcs_bai_path\na\tgs://b/a.bam\tgs://b/a.bam.bai\n")
    with pytest.raises(SystemExit) as exc:
        mpt.read_bams_table(str(p), "gcs_bam_path", "gcs_bai_path")
    assert "entity:<type>_id" in str(exc.value)


def test_reader_collects_paths_in_file_order_and_flags_duplicates(tmp_path: Path):
    p = tmp_path / "t.txt"
    p.write_text("a\tgs://b/early.vcf.gz\na\tgs://b/late.vcf.gz\nb\tgs://b/x.vcf.gz\n")
    paths, dups = mpt.read_trgt_table(str(p))
    assert paths == {"a": ["gs://b/early.vcf.gz", "gs://b/late.vcf.gz"], "b": ["gs://b/x.vcf.gz"]}
    assert dups == ["a"]


def test_repeated_identical_path_is_not_a_duplicate(tmp_path: Path):
    p = tmp_path / "t.txt"
    p.write_text("a\tgs://b/x.vcf.gz\na\tgs://b/x.vcf.gz\n")
    paths, dups = mpt.read_trgt_table(str(p))
    assert paths == {"a": ["gs://b/x.vcf.gz"]} and dups == []


def test_resolve_picks_first_or_last(tmp_path: Path):
    paths = {"a": ["early", "late"]}
    assert mpt.resolve(paths, "first") == {"a": "early"}
    assert mpt.resolve(paths, "last") == {"a": "late"}
    assert mpt.resolve(paths, "fail") == {"a": "early"}, "fail resolves like first; run() guards"
