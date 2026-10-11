#!/usr/bin/env python3
"""Tests for felix/scripts/compare_ru_arms.py (four-arm comparison table)."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from compare_ru_arms import lambda_from_p, load, main  # noqa: E402

NUM_ANCS = 3


def _results(tmp: Path, name: str, pooled: list[float], cct: list[float], n: int) -> Path:
    """FELIX admixed single-variant output: keys, per-ancestry block, joint Ps."""
    df = pd.DataFrame({
        "CHR": ["chr22"] * n,
        "POS": list(range(1000, 1000 + n)),
        "MarkerID": [f"TR{i}" for i in range(n)],
        "Allele1": ["N"] * n,
        "Allele2": ["<TR>"] * n,
    })
    for k in range(1, NUM_ANCS + 1):
        df[f"AC_Allele2_anc{k}"] = 30 + k
        df[f"BETA_anc{k}"] = 0.1 * k
        df[f"p.value_anc{k}"] = 0.5
        df[f"p.value_c_anc{k}"] = [min(1.0, 0.5 * (i + 1) / n) for i in range(n)]
    df["AC_Allele2_ancALL"] = list(range(20, 20 + n))
    df["BETA_ancALL"] = 0.2
    df["p.value_ancALL"] = pooled
    df["P_het_admixed_c"] = 0.4
    df["P_hom_admixed_c"] = pooled
    df["P_cct_admixed_c"] = cct
    path = tmp / f"{name}.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_picks_felix_columns_and_carrier_counts():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        p = _results(tmp, "bia", [0.5] * 4, [0.5] * 4, 4)
        out, prov = load(p, "biallelic", NUM_ANCS)
    assert prov["pooled_p_column"] == "p.value_ancALL"
    assert prov["la_p_column"] == "P_cct_admixed_c"
    assert prov["beta_column"] == "BETA_ancALL"
    # Conditional per-ancestry columns win over the unconditional ones.
    assert prov["per_ancestry_p_columns"] == {
        "anc1": "p.value_c_anc1", "anc2": "p.value_c_anc2", "anc3": "p.value_c_anc3"
    }
    assert list(out["n_carrier_hap"]) == [20, 21, 22, 23]
    assert "P_minanc_biallelic" in out.columns


def test_contrasts_count_gained_and_lost():
    n = 50
    base = [0.5] * n
    # Locus 0: significant only with dosage. Locus 1: only biallelic (a loss).
    bia_pooled, dos_pooled = list(base), list(base)
    bia_cct, dos_cct = list(base), list(base)
    dos_pooled[0] = 1e-12
    bia_pooled[1] = 1e-12
    dos_cct[0] = 1e-20
    dos_cct[2] = 1e-20
    bia_cct[1] = 1e-20
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        bia = _results(tmp, "bia", bia_pooled, bia_cct, n)
        dos = _results(tmp, "dos", dos_pooled, dos_cct, n)
        argv = [
            "compare_ru_arms.py",
            "--dosage-tsv", str(dos), "--biallelic-tsv", str(bia),
            "--phenotype", "T2D", "--num-ancs", str(NUM_ANCS),
            "--out-tsv", str(tmp / "arms.tsv"),
            "--out-summary-json", str(tmp / "arms.json"),
            "--out-md", str(tmp / "arms.md"),
        ]
        old, sys.argv = sys.argv, argv
        try:
            assert main() == 0
        finally:
            sys.argv = old
        summary = json.loads((tmp / "arms.json").read_text())
        table = pd.read_csv(tmp / "arms.tsv", sep="\t")
        md = (tmp / "arms.md").read_text()

    assert summary["n_loci_shared"] == n
    assert summary["n_loci_biallelic_only"] == 0 and summary["n_loci_dosage_only"] == 0
    assert summary["arms"]["pooled_biallelic"]["n_significant"] == 1
    assert summary["arms"]["pooled_dosage"]["n_significant"] == 1
    assert summary["arms"]["la_dosage"]["n_significant"] == 2
    coding = summary["contrasts"]["pooled_biallelic->pooled_dosage"]
    assert (coding["gained"], coding["lost"], coding["both"]) == (1, 1, 0)
    ancestry = summary["contrasts"]["pooled_dosage->la_dosage"]
    assert (ancestry["gained"], ancestry["lost"], ancestry["both"]) == (1, 0, 1)
    # Table is sorted by the strongest arm and keeps all four p-columns.
    assert list(table.columns[:4]) == ["CHR", "POS", "MarkerID", "n_carrier_hap"]
    assert table["MarkerID"].iloc[0] in {"TR0", "TR2"}
    assert {"P_pooled_biallelic", "P_la_biallelic", "P_pooled_dosage", "P_la_dosage"} <= set(table.columns)
    assert "not 1-df" in md and "| la_dosage |" in md


def test_unmatched_marker_sets_are_reported_not_silent():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        bia = _results(tmp, "bia", [0.5] * 12, [0.5] * 12, 12)
        dos = _results(tmp, "dos", [0.5] * 10, [0.5] * 10, 10)
        out_bia, _ = load(bia, "biallelic", NUM_ANCS)
        out_dos, _ = load(dos, "dosage", NUM_ANCS)
        merged = out_bia.merge(out_dos, on=["CHR", "POS", "MarkerID"], how="inner")
        assert len(merged) == 10 and len(out_bia) - len(merged) == 2


def test_missing_p_column_fails_loudly():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        p = _results(tmp, "bia", [0.5] * 4, [0.5] * 4, 4)
        df = pd.read_csv(p, sep="\t").drop(columns=["P_cct_admixed_c"])
        df.to_csv(p, sep="\t", index=False)
        try:
            load(p, "biallelic", NUM_ANCS)
        except SystemExit as exc:
            assert "local-ancestry p-value" in str(exc)
        else:
            raise AssertionError("expected missing-column failure")


def test_lambda_matches_a_known_null():
    rng = np.random.default_rng(0)
    lam = lambda_from_p(rng.uniform(size=20000))
    assert abs(lam - 1.0) < 0.05, lam
    # Inflated input must push lambda well above 1.
    inflated = rng.uniform(size=20000) ** 2
    assert lambda_from_p(inflated) > 1.3


if __name__ == "__main__":
    test_picks_felix_columns_and_carrier_counts()
    test_contrasts_count_gained_and_lost()
    test_unmatched_marker_sets_are_reported_not_silent()
    test_missing_p_column_fails_loudly()
    test_lambda_matches_a_known_null()
    print("ok")
