"""Tests for estimate_repeat_exposure.py."""

from __future__ import annotations

import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("ere", _HERE / "estimate_repeat_exposure.py")
ere = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ere)

HEADER = "SAMPLE\teas\tamr\teur\tafr\tsas\n"


def write_anc(tmp_path: Path, rows: list[str], gz: bool = True) -> Path:
    p = tmp_path / ("global.anc.gz" if gz else "global.anc")
    data = HEADER + "".join(rows)
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(data)
    else:
        p.write_text(data)
    return p


def run_script(anc: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(_HERE / "estimate_repeat_exposure.py"), "--global-anc", str(anc), *extra],
        capture_output=True, text=True)


def test_single_ancestry_has_no_exposure():
    """A participant with one ancestry can never have discordant haplotypes."""
    assert ere.discordance(np.array([[0.0, 0.0, 0.0, 1.0, 0.0]]))[0] == pytest.approx(0.0)


def test_even_two_way_admixture():
    """50/50 gives 1 - (0.25 + 0.25) = 0.5."""
    assert ere.discordance(np.array([[0.0, 0.0, 0.5, 0.5, 0.0]]))[0] == pytest.approx(0.5)


def test_even_five_way_admixture():
    """Equal shares over k populations give 1 - 1/k."""
    p = np.full((1, 5), 0.2)
    assert ere.discordance(p)[0] == pytest.approx(0.8)


def test_proportions_are_renormalised():
    """FLARE rows need not sum to exactly 1; the measure must not care."""
    a = ere.discordance(np.array([[0.0, 0.0, 0.5, 0.5, 0.0]]))[0]
    b = ere.discordance(np.array([[0.0, 0.0, 0.25, 0.25, 0.0]]))[0]
    assert a == pytest.approx(b)


def test_all_zero_row_is_not_a_division_error():
    assert ere.discordance(np.zeros((1, 5)))[0] == pytest.approx(1.0)


def test_discordance_is_bounded():
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(5), size=500)
    d = ere.discordance(p)
    assert (d >= 0).all() and (d <= 1 - 1 / 5 + 1e-12).all()


def test_header_sample_column_is_renamed(tmp_path: Path):
    anc = write_anc(tmp_path, ["s1\t0\t0\t0\t1\t0\n"])
    df, pops = ere.read_global_anc(str(anc))
    assert list(df.columns)[0] == "sample"
    assert pops == ["eas", "amr", "eur", "afr", "sas"]
    assert df["sample"].iloc[0] == "s1"


def test_sample_ids_stay_strings(tmp_path: Path):
    """Person ids are numeric-looking; they must not become ints or floats."""
    anc = write_anc(tmp_path, ["1000612\t0\t0\t0\t0.999\t0\n"])
    df, _ = ere.read_global_anc(str(anc))
    assert df["sample"].iloc[0] == "1000612"


def test_end_to_end_reports_expected_counts(tmp_path: Path):
    anc = write_anc(tmp_path, [
        "s_afr\t0\t0\t0\t1\t0\n",          # no exposure
        "s_mix\t0\t0\t0.5\t0.5\t0\n",      # 0.5
    ])
    out = tmp_path / "per_sample.tsv"
    proc = run_script(anc, "--length-het-loci", "2431", "--window-loci", "14618",
                      "--catalog-loci", "4439813", "--out", str(out))
    assert proc.returncode == 0, proc.stderr
    assert "participants:              2" in proc.stdout
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[-2:] == ["p_anc_discordant", "exposed_loci"]
    vals = {l.split("\t")[0]: float(l.split("\t")[-1]) for l in lines[1:]}
    assert vals["s_afr"] == pytest.approx(0.0)
    # 0.5 * (2431/14618) * 4,439,813
    assert vals["s_mix"] == pytest.approx(0.5 * (2431 / 14618) * 4439813, rel=1e-6)


def test_uncompressed_input_also_works(tmp_path: Path):
    anc = write_anc(tmp_path, ["s1\t0\t0\t0.5\t0.5\t0\n"], gz=False)
    assert run_script(anc).returncode == 0
