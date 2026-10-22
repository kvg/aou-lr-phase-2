"""Tests for make_trgt_split_table.py."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("mtst", _HERE / "make_trgt_split_table.py")
mtst = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(mtst)

LAI_HEADER = "entity:flare_lai_prod_id\tchrom\tgt_vcf\tgt_vcf_index\n"


def lai_table(tmp_path: Path, chroms=("chr1", "chr2")) -> Path:
    p = tmp_path / "lai.tsv"
    rows = [LAI_HEADER]
    for c in chroms:
        rows.append(f"aou.{c}\t{c}\tgs://b/{c}.vcf.gz\tgs://b/{c}.vcf.gz.tbi\n")
    p.write_text("".join(rows))
    return p


def sample_list(tmp_path: Path, n: int, header: bool = False) -> Path:
    p = tmp_path / "samples.tsv"
    rows = ["entity:trgt_phase_id\tsample\n"] if header else []
    rows += [f"100{i}\t100{i}\n" if header else f"100{i}\n" for i in range(n)]
    p.write_text("".join(rows))
    return p


def run_script(tmp_path: Path, lai: Path, samples: Path, *extra: str):
    table, batches = tmp_path / "table.tsv", tmp_path / "batches.tsv"
    proc = subprocess.run(
        [sys.executable, str(_HERE / "make_trgt_split_table.py"),
         "--lai-table", str(lai), "--samples", str(samples),
         "--out-base", "gs://out/split", "--out-table", str(table),
         "--out-batches", str(batches), *extra],
        capture_output=True, text=True)
    return proc, table, batches


# --- batching ---------------------------------------------------------------

def test_batches_cover_every_sample_exactly_once():
    samples = [f"s{i}" for i in range(1050)]
    rows = mtst.batch_rows(samples, 500)
    seen = [s for _, joined in rows for s in joined.split(",")]
    assert seen == samples, "order preserved, nothing dropped or duplicated"
    assert len(rows) == 3
    assert [len(j.split(",")) for _, j in rows] == [500, 500, 50]


def test_batch_names_are_zero_padded_consistently():
    """batch-0 and batch-10 sorting apart would scramble the scatter's file order."""
    rows = mtst.batch_rows([f"s{i}" for i in range(2500)], 500)
    names = [n for n, _ in rows]
    assert names == sorted(names), "lexical order must match numeric order"
    assert names[0] == "batch-0"


def test_single_batch_when_cohort_fits():
    rows = mtst.batch_rows(["a", "b"], 500)
    assert rows == [("batch-0", "a,b")]


# --- inputs -----------------------------------------------------------------

def test_terra_export_header_is_skipped(tmp_path: Path):
    assert mtst.read_samples(str(sample_list(tmp_path, 3, header=True))) == ["1000", "1001", "1002"]


def test_bare_id_list_is_accepted(tmp_path: Path):
    assert mtst.read_samples(str(sample_list(tmp_path, 3))) == ["1000", "1001", "1002"]


def test_duplicate_sample_ids_are_collapsed(tmp_path: Path):
    p = tmp_path / "s.txt"
    p.write_text("1000\n1001\n1000\n")
    assert mtst.read_samples(str(p)) == ["1000", "1001"]


def test_missing_lai_column_names_the_columns(tmp_path: Path):
    p = tmp_path / "lai.tsv"
    p.write_text("entity:x_id\tchrom\tsomething_else\naou.chr1\tchr1\tv\n")
    with pytest.raises(SystemExit, match="gt_vcf"):
        mtst.read_lai_table(str(p), "chrom", "gt_vcf", "gt_vcf_index")


def test_rows_with_an_empty_vcf_cell_are_skipped(tmp_path: Path):
    p = tmp_path / "lai.tsv"
    p.write_text(LAI_HEADER + "aou.chr1\tchr1\t\t\naou.chr2\tchr2\tgs://b/2\tgs://b/2.tbi\n")
    got = mtst.read_lai_table(str(p), "chrom", "gt_vcf", "gt_vcf_index")
    assert set(got) == {"chr2"}, "an unfilled FLARE row is not a usable split input"


# --- end to end -------------------------------------------------------------

def test_table_has_one_row_per_chromosome(tmp_path: Path):
    """NOT sub-chromosomal: a shard writes per-sample files for its region, so
    sharding would multiply the file count by the shard count."""
    proc, table, batches = run_script(
        tmp_path, lai_table(tmp_path, ("chr1", "chr2", "chr3")), sample_list(tmp_path, 10))
    assert proc.returncode == 0, proc.stderr
    lines = table.read_text().splitlines()
    assert lines[0].split("\t") == [
        "entity:trgt_split_regions_id", "region", "joint_vcf", "joint_vcf_idx",
        "gcs_output_dir", "output_tag", "sample_batches_tsv"]
    rows = [l.split("\t") for l in lines[1:]]
    assert [r[0] for r in rows] == ["chr1", "chr2", "chr3"]
    assert [r[1] for r in rows] == ["chr1", "chr2", "chr3"], "region is the whole contig"


def test_output_tag_makes_paths_deterministic(tmp_path: Path):
    proc, table, _ = run_script(tmp_path, lai_table(tmp_path, ("chr7",)),
                                sample_list(tmp_path, 3), "--output-prefix", "phased")
    assert proc.returncode == 0, proc.stderr
    row = dict(zip(table.read_text().splitlines()[0].split("\t"),
                   table.read_text().splitlines()[1].split("\t")))
    assert row["gcs_output_dir"] == "gs://out/split/chr7"
    assert row["output_tag"] == "phased.chr7"
    assert "vcf.gz" in proc.stdout and ".bcf" not in proc.stdout


def test_batches_file_uses_the_documented_three_column_form(tmp_path: Path):
    proc, _, batches = run_script(tmp_path, lai_table(tmp_path), sample_list(tmp_path, 7),
                                  "--batch-size", "3")
    assert proc.returncode == 0, proc.stderr
    rows = [l.split("\t") for l in batches.read_text().splitlines()]
    assert len(rows) == 3
    assert all(len(r) == 3 and r[1] == "-" for r in rows), "samples, '-', batch name"
    assert rows[0][2] == "batch-0"
    assert [len(r[0].split(",")) for r in rows] == [3, 3, 1]


def test_autosomes_absent_from_the_lai_table_are_reported(tmp_path: Path):
    proc, table, _ = run_script(tmp_path, lai_table(tmp_path, ("chr1",)), sample_list(tmp_path, 2))
    assert proc.returncode == 0, proc.stderr
    assert len(table.read_text().splitlines()) == 2
    assert "autosomes absent" in proc.stderr
    assert "chr22" in proc.stderr


def test_non_autosomes_are_ignored(tmp_path: Path):
    proc, table, _ = run_script(tmp_path, lai_table(tmp_path, ("chr1", "chrX", "chrM")),
                                sample_list(tmp_path, 2))
    assert proc.returncode == 0, proc.stderr
    assert [l.split("\t")[0] for l in table.read_text().splitlines()[1:]] == ["chr1"]
    assert "chrM" in proc.stderr and "chrX" in proc.stderr


def test_local_batches_path_warns(tmp_path: Path):
    """Terra needs a gs:// path; a local one would resolve to nothing."""
    proc, _, _ = run_script(tmp_path, lai_table(tmp_path), sample_list(tmp_path, 2))
    assert "Terra needs a gs:// path" in proc.stderr


def test_batches_uri_is_used_in_the_table(tmp_path: Path):
    proc, table, _ = run_script(tmp_path, lai_table(tmp_path), sample_list(tmp_path, 2),
                                "--batches-uri", "gs://b/batches.tsv")
    assert proc.returncode == 0, proc.stderr
    assert all(l.split("\t")[-1] == "gs://b/batches.tsv"
               for l in table.read_text().splitlines()[1:])
    assert "Terra needs a gs:// path" not in proc.stderr


def test_empty_sample_list_is_refused(tmp_path: Path):
    p = tmp_path / "s.txt"
    p.write_text("")
    proc, _, _ = run_script(tmp_path, lai_table(tmp_path), p)
    assert proc.returncode != 0
    assert "no sample ids" in proc.stderr
