"""Tests for make_trgt_phase_table.py."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("mtpt", _HERE / "make_trgt_phase_table.py")
mtpt = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(mtpt)

BI = "gs://b/CCS/terra-d45cc5ca/outputs/{}_trgt.vcf.gz"
UW = "gs://b/CCS/terra-f6671367/outputs/{}_trgt.vcf.gz"


def write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def bams_table(tmp_path: Path, ids: list[str], blank: set[str] = frozenset()) -> Path:
    lines = ["entity:sample_id\tgcs_bam_path\tgcs_bai_path"]
    for i in ids:
        bam = "" if i in blank else f"gs://x/{i}.bam"
        bai = "" if i in blank else f"gs://x/{i}.bam.bai"
        lines.append(f"{i}\t{bam}\t{bai}")
    return write(tmp_path / "bams.tsv", "\n".join(lines) + "\n")


def run_script(tmp_path: Path, bams: Path, trgt: Path, samples: Path,
               *extra: str) -> tuple[subprocess.CompletedProcess, Path]:
    out = tmp_path / "trgt_phase.tsv"
    proc = subprocess.run(
        [sys.executable, str(_HERE / "make_trgt_phase_table.py"),
         "--bams-table", str(bams), "--trgt-table", str(trgt),
         "--sample-list", str(samples), "--bam-column", "gcs_bam_path",
         "--bai-column", "gcs_bai_path", "--out", str(out), *extra],
        capture_output=True, text=True)
    return proc, out


# --- site mapping -----------------------------------------------------------

def test_submission_maps_to_site_by_majority_vote():
    rows = [("p1", BI.format("p1")), ("p2", BI.format("p2")), ("p3", UW.format("p3"))]
    sites = {"p1": ["BI"], "p2": ["BI"], "p3": ["UW"]}
    mapping, notes = mtpt.submission_sites(rows, sites, 0.95)
    assert mapping == {"terra-d45cc5ca": "BI", "terra-f6671367": "UW"}
    assert any("terra-d45cc5ca" in n and "BI" in n for n in notes)


def test_ambiguous_people_do_not_vote():
    """A person at two sites is what the mapping disambiguates; it cannot also vote."""
    rows = [("dup", BI.format("dup")), ("p1", BI.format("p1"))]
    sites = {"dup": ["BI", "UW"], "p1": ["BI"]}
    _, notes = mtpt.submission_sites(rows, sites, 0.95)
    assert "1/1" in notes[0], notes


def test_impure_submission_is_refused():
    rows = [("p1", BI.format("p1")), ("p2", BI.format("p2")), ("p3", BI.format("p3"))]
    sites = {"p1": ["BI"], "p2": ["UW"], "p3": ["UW"]}
    with pytest.raises(SystemExit, match="does not map cleanly"):
        mtpt.submission_sites(rows, sites, 0.95)


# --- keying -----------------------------------------------------------------

def test_single_site_person_ignores_the_submission():
    """An unambiguous person takes their own site even if the batch says otherwise."""
    rows = [("p1", UW.format("p1"))]
    trgt, tally = mtpt.key_by_sample(rows, {"p1": ["BI"]}, {"terra-f6671367": "UW"})
    assert trgt == {"p1_BI": ("p1", UW.format("p1"))}
    assert not tally


def test_cross_site_person_becomes_two_rows():
    """The 28 repeats are two samples, not a collision."""
    rows = [("dup", BI.format("dup")), ("dup", UW.format("dup"))]
    sites = {"dup": ["BI", "UW"]}
    mapping = {"terra-d45cc5ca": "BI", "terra-f6671367": "UW"}
    trgt, tally = mtpt.key_by_sample(rows, sites, mapping)
    assert set(trgt) == {"dup_BI", "dup_UW"}
    assert trgt["dup_BI"][1] == BI.format("dup")
    assert trgt["dup_UW"][1] == UW.format("dup")
    assert not tally


def test_person_absent_from_sample_list_is_tallied():
    trgt, tally = mtpt.key_by_sample([("ghost", BI.format("ghost"))], {}, {})
    assert trgt == {}
    assert tally["person_not_in_sample_list"] == 1


def test_submission_site_not_among_persons_sites_is_tallied():
    rows = [("dup", "gs://b/CCS/terra-aaaaaaaa/outputs/dup_trgt.vcf.gz")]
    trgt, tally = mtpt.key_by_sample(rows, {"dup": ["BI", "UW"]}, {"terra-aaaaaaaa": "HA"})
    assert trgt == {}
    assert tally["submission_site_not_a_site_of_person"] == 1


# --- sample list parsing ----------------------------------------------------

def test_membership_export_is_accepted(tmp_path: Path):
    p = write(tmp_path / "m.tsv",
              "membership:sample_set_id\tsample\nset1\tp1_BI\nset1\tp2_UW\n")
    assert mtpt.read_sample_list(str(p)) == {"p1": ["BI"], "p2": ["UW"]}


def test_bare_id_list_is_accepted(tmp_path: Path):
    p = write(tmp_path / "ids.txt", "p1_BI\np1_UW\np2_HA\n")
    assert mtpt.read_sample_list(str(p)) == {"p1": ["BI", "UW"], "p2": ["HA"]}


def test_sample_list_without_site_suffixes_is_refused(tmp_path: Path):
    p = write(tmp_path / "ids.txt", "p1\np2\n")
    with pytest.raises(SystemExit, match="no '<person>_<SITE>' ids"):
        mtpt.read_sample_list(str(p))


# --- end to end -------------------------------------------------------------

def test_entity_id_is_the_sample_and_sample_column_is_the_person(tmp_path: Path):
    """bcftools -s must find the PERSON in the joint callset, not the sample id."""
    trgt = write(tmp_path / "t.txt", f"p1\t{BI.format('p1')}\n")
    samples = write(tmp_path / "s.txt", "p1_BI\n")
    proc, out = run_script(tmp_path, bams_table(tmp_path, ["p1_BI"]), trgt, samples)
    assert proc.returncode == 0, proc.stderr
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == [
        "entity:trgt_phase_id", "sample", "bam", "bam_index", "trgt_vcf", "trgt_vcf_index"]
    row = lines[1].split("\t")
    assert row[0] == "p1_BI", "entity id is the sample"
    assert row[1] == "p1", "sample column is the bare person id"
    assert row[5] == row[4] + ".tbi"


def test_only_the_intersection_is_written(tmp_path: Path):
    trgt = write(tmp_path / "t.txt",
                 f"p1\t{BI.format('p1')}\nq2\t{BI.format('q2')}\n")
    samples = write(tmp_path / "s.txt", "p1_BI\nq1_BI\nq2_BI\n")
    bams = bams_table(tmp_path, ["p1_BI", "q1_BI"])
    proc, out = run_script(tmp_path, bams, trgt, samples)
    assert proc.returncode == 0, proc.stderr
    assert [l.split("\t")[0] for l in out.read_text().splitlines()[1:]] == ["p1_BI"]
    assert "BAM but no TRGT call:         1" in proc.stdout
    assert "TRGT call but no usable BAM:  1" in proc.stdout


def test_replicate_pair_survives_end_to_end(tmp_path: Path):
    trgt = write(tmp_path / "t.txt",
                 f"dup\t{BI.format('dup')}\ndup\t{UW.format('dup')}\n"
                 f"p1\t{BI.format('p1')}\n")
    samples = write(tmp_path / "s.txt", "dup_BI\ndup_UW\np1_BI\n")
    bams = bams_table(tmp_path, ["dup_BI", "dup_UW", "p1_BI"])
    proc, out = run_script(tmp_path, bams, trgt, samples)
    assert proc.returncode == 0, proc.stderr
    rows = {l.split("\t")[0]: l.split("\t") for l in out.read_text().splitlines()[1:]}
    assert set(rows) == {"dup_BI", "dup_UW", "p1_BI"}
    assert rows["dup_BI"][1] == rows["dup_UW"][1] == "dup", "same person, two samples"
    assert rows["dup_BI"][4] != rows["dup_UW"][4], "different TRGT calls"
    assert "with both site samples kept:  1" in proc.stdout


def test_blank_bam_cell_excludes_the_sample(tmp_path: Path):
    trgt = write(tmp_path / "t.txt", f"p1\t{BI.format('p1')}\n")
    samples = write(tmp_path / "s.txt", "p1_BI\n")
    bams = bams_table(tmp_path, ["p1_BI"], blank={"p1_BI"})
    proc, out = run_script(tmp_path, bams, trgt, samples)
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().splitlines()[1:] == []
    assert "missing a BAM or index cell:  1" in proc.stdout


def test_rows_follow_the_bam_table_order(tmp_path: Path):
    trgt = write(tmp_path / "t.txt",
                 "".join(f"p{i}\t{BI.format(f'p{i}')}\n" for i in (1, 2, 3)))
    samples = write(tmp_path / "s.txt", "p1_BI\np2_BI\np3_BI\n")
    bams = bams_table(tmp_path, ["p3_BI", "p1_BI", "p2_BI"])
    proc, out = run_script(tmp_path, bams, trgt, samples)
    assert proc.returncode == 0, proc.stderr
    assert [l.split("\t")[0] for l in out.read_text().splitlines()[1:]] == [
        "p3_BI", "p1_BI", "p2_BI"]


def test_non_terra_bams_table_is_refused(tmp_path: Path):
    bams = write(tmp_path / "bams.tsv", "sample\tgcs_bam_path\tgcs_bai_path\np1_BI\ta\tb\n")
    trgt = write(tmp_path / "t.txt", f"p1\t{BI.format('p1')}\n")
    samples = write(tmp_path / "s.txt", "p1_BI\n")
    proc, _ = run_script(tmp_path, bams, trgt, samples)
    assert proc.returncode != 0
    assert "expected 'entity:<type>_id'" in proc.stderr


def test_malformed_trgt_row_is_refused(tmp_path: Path):
    trgt = write(tmp_path / "t.txt", "p1\n")
    samples = write(tmp_path / "s.txt", "p1_BI\n")
    proc, _ = run_script(tmp_path, bams_table(tmp_path, ["p1_BI"]), trgt, samples)
    assert proc.returncode != 0
    assert "expected 'id<TAB>uri'" in proc.stderr


def test_repeated_identical_path_is_collapsed(tmp_path: Path):
    """The same row twice is one sample, not a replicate pair."""
    rows = mtpt.read_trgt_table(str(write(
        tmp_path / "t.txt", f"p1\t{BI.format('p1')}\np1\t{BI.format('p1')}\n")))
    assert rows == [("p1", BI.format("p1"))]


def test_person_id_containing_an_underscore_parses(tmp_path: Path):
    """The site is a suffix, so the split must come from the right."""
    p = write(tmp_path / "ids.txt", "odd_name_BI\n")
    assert mtpt.read_sample_list(str(p)) == {"odd_name": ["BI"]}
