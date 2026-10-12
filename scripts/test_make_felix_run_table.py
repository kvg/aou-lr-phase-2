"""Tests for scripts/make_felix_run_table.py."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "make_felix_run_table.py"
TEMPLATE = ROOT / "felix" / "configs" / "terra" / "felix_runs.entity.tsv"
ID_COL = "entity:flare_lai_prod_set_id"
BUCKET = "gs://fc-secure-test"


def run(tmp_path, names, *extra, phenotypes=None):
    ph = phenotypes or tmp_path / "phenos.txt"
    if phenotypes is None:
        ph.write_text("\n".join(names) + "\n")
    out = tmp_path / "out"
    cp = subprocess.run(
        [sys.executable, str(SCRIPT), "--phenotypes", str(ph), "--template", str(TEMPLATE),
         "--template-row", "felix.wgs.limited.v1", "--bucket", BUCKET,
         "--out-dir", str(out), *extra],
        capture_output=True, text=True)
    return cp, out


def read(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_batches_cover_every_phenotype_exactly_once(tmp_path):
    names = [f"PH_{i:03d}" for i in range(943)]
    cp, out = run(tmp_path, names, "--batch-size", "95")
    assert cp.returncode == 0, cp.stderr
    files = sorted((out / "phenotypes").glob("*.txt"))
    assert len(files) == 10
    seen = [ln for f in files for ln in f.read_text().split()]
    assert seen == names                      # order preserved, nothing lost or repeated


def test_rows_differ_only_in_id_phenotypes_and_description(tmp_path):
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(20)], "--batch-size", "7")
    assert cp.returncode == 0, cp.stderr
    rows = read(out / "felix_runs.entity.tsv")
    assert len(rows) == 3
    varying = {c for c in rows[0] if len({r[c] for r in rows}) > 1}
    assert varying == {ID_COL, "selected_phenotypes", "run_description"}
    template = {r[ID_COL]: r for r in read(TEMPLATE)}["felix.wgs.limited.v1"]
    for r in rows:
        for c in r:
            if c not in varying:
                assert r[c] == template[c], c


def test_phenotype_cells_point_at_the_written_files(tmp_path):
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(5)], "--batch-size", "2",
                  "--phenotype-prefix", "felix/ph")
    assert cp.returncode == 0, cp.stderr
    for r in read(out / "felix_runs.entity.tsv"):
        uri = r["selected_phenotypes"]
        assert uri.startswith(f"{BUCKET}/felix/ph/")
        assert (out / "phenotypes" / uri.rsplit("/", 1)[1]).is_file()


def test_membership_is_two_columns_and_chromosome_ordered(tmp_path):
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(5)], "--batch-size", "2")
    assert cp.returncode == 0, cp.stderr
    lines = (out / "felix_runs.membership.tsv").read_text().rstrip("\n").split("\n")
    assert lines[0] == "membership:flare_lai_prod_set_id\tflare_lai_prod"
    body = [ln.split("\t") for ln in lines[1:]]
    assert all(len(f) == 2 for f in body)
    assert len(body) == 3 * 22
    first = [f[1].split(".")[-1] for f in body[:22]]
    assert first == [f"chr{i}" for i in range(1, 23)]


def test_pheno_cov_override_applies_to_every_row(tmp_path):
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(5)], "--batch-size", "2",
                  "--pheno-cov", "gs://b/new_pheno_cov.tsv")
    assert cp.returncode == 0, cp.stderr
    assert {r["pheno_cov"] for r in read(out / "felix_runs.entity.tsv")} == {
        "gs://b/new_pheno_cov.tsv"}


def test_batches_are_balanced_not_sliced(tmp_path):
    """--batch-size is a ceiling. Slicing at it would leave a batch of two."""
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(10)], "--batch-size", "4")
    assert cp.returncode == 0, cp.stderr
    sizes = [len(f.read_text().split()) for f in sorted((out / "phenotypes").glob("*.txt"))]
    assert sizes == [4, 3, 3]


def test_no_batch_of_one_at_the_real_scale(tmp_path):
    """951 phenotypes at 95 slices to a final batch of a single phenotype."""
    cp, out = run(tmp_path, [f"PH_{i:03d}" for i in range(951)], "--batch-size", "95")
    assert cp.returncode == 0, cp.stderr
    sizes = [len(f.read_text().split()) for f in sorted((out / "phenotypes").glob("*.txt"))]
    assert sum(sizes) == 951
    assert max(sizes) <= 95 and max(sizes) - min(sizes) <= 1, sizes
    assert sizes == [87] * 5 + [86] * 6, sizes      # 951 = 5*87 + 6*86


@pytest.mark.parametrize("n,size", [(1, 95), (95, 95), (96, 95), (943, 95), (951, 95), (650, 95)])
def test_batching_is_exhaustive_ordered_and_capped(tmp_path, n, size):
    names = [f"PH_{i:04d}" for i in range(n)]
    cp, out = run(tmp_path, names, "--batch-size", str(size))
    assert cp.returncode == 0, cp.stderr
    files = sorted((out / "phenotypes").glob("*.txt"))
    seen = [x for f in files for x in f.read_text().split()]
    assert seen == names
    sizes = [len(f.read_text().split()) for f in files]
    assert max(sizes) <= size and min(sizes) >= 1
    assert len(files) == -(-n // size)


@pytest.mark.parametrize("names,extra,msg", [
    ([], [], "empty"),
    (["A", "B", "A"], [], "duplicate"),
    (["A"], ["--batch-size", "0"], "batch-size"),
])
def test_bad_input_is_refused(tmp_path, names, extra, msg):
    cp, _ = run(tmp_path, names, *extra)
    assert cp.returncode != 0
    assert msg in cp.stderr


def test_auto_created_row_is_refused(tmp_path):
    ph = tmp_path / "p.txt"
    ph.write_text("A\n")
    cp = subprocess.run(
        [sys.executable, str(SCRIPT), "--phenotypes", str(ph), "--template", str(TEMPLATE),
         "--template-row", "PropagateFlareAncestry_2026-10-10T22-08-23",
         "--bucket", BUCKET, "--out-dir", str(tmp_path / "o")],
        capture_output=True, text=True)
    assert cp.returncode != 0
    assert "not in" in cp.stderr or "no values" in cp.stderr


def test_set_overrides_any_column(tmp_path):
    cp, out = run(tmp_path, [f"PH_{i}" for i in range(5)], "--batch-size", "2",
                  "--set", "min_mac=20",
                  "--set", "covariate_columns=gs://b/covariate_columns_full.txt")
    assert cp.returncode == 0, cp.stderr
    rows = read(out / "felix_runs.entity.tsv")
    assert {r["min_mac"] for r in rows} == {"20"}
    assert {r["covariate_columns"] for r in rows} == {"gs://b/covariate_columns_full.txt"}


@pytest.mark.parametrize("bad,msg", [
    ("not_a_column=1", "not a column"),
    ("min_mac", "COLUMN=VALUE"),
    ("entity:flare_lai_prod_set_id=x", "cannot change"),
])
def test_bad_set_is_refused(tmp_path, bad, msg):
    cp, _ = run(tmp_path, ["A", "B"], "--batch-size", "1", "--set", bad)
    assert cp.returncode != 0 and msg in cp.stderr
