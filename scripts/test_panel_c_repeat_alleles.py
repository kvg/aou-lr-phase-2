#!/usr/bin/env python3
"""Tests for scripts/panel_c_repeat_alleles.py."""

from __future__ import annotations

import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aggregate_repeat_loci import ALLELE_COLUMNS  # noqa: E402
from panel_c_repeat_alleles import expected_alleles, log_factorials, main  # noqa: E402


def test_expected_alleles_matches_enumeration():
    counts = [np.array([3, 2, 1]), np.array([6])]
    haps = [0, 0, 0, 1, 1, 2]
    for n in (1, 2, 3, 5):
        brute = np.mean([len(set(c)) for c in itertools.combinations(haps, n)])
        got = expected_alleles(counts, np.array([n]), log_factorials(6))
        assert got[0, 0] == pytest.approx(brute)
        assert got[1, 0] == pytest.approx(1.0)


def _table(path: Path, rows: list[tuple[str, int, int, list[tuple[int, int]]]]) -> str:
    with path.open("w") as fh:
        fh.write("\t".join(ALLELE_COLUMNS) + "\n")
        for lid, n_hap, n_missing, alleles in rows:
            for bp, c in alleles:
                fh.write(f"{lid}\tchr1\t100\t130\t3\tCAG\t1\t2\t{n_hap}\t{n_missing}\t{bp}\t{bp / 3:g}\t{c}\n")
    return str(path)


def test_rarefaction_examples_and_new_alleles(tmp_path: Path):
    p1 = _table(tmp_path / "p1.tsv", [
        ("A", 20, 0, [(0, 10), (3, 6), (6, 4)]),
        ("B", 20, 0, [(0, 19), (-3, 1)]),
        ("LOW", 10, 10, [(0, 10)]),  # below the call-rate floor
    ])
    p2 = _table(tmp_path / "p2.tsv", [
        ("A", 200, 0, [(-3, 10), (0, 100), (3, 50), (6, 30), (9, 8), (12, 2)]),
        ("C", 198, 2, [(0, 190), (30, 8)]),
    ])
    plvi = tmp_path / "plvi.tsv"
    plvi.write_text(
        "trid\tchrom\tstart\tend\tmotif\tmotif_len\tmotif_group\tn_hap\tlps_mean\tplvi\tplvi_pct_in_group\tn_in_group\n"
        "C\tchr1\t0\t1\tCAG\t3\t3\t200\t10\t9\t100\t3\n"
        "A\tchr1\t0\t1\tCAG\t3\t3\t200\t10\t5\t66.7\t3\n"
        "B\tchr1\t0\t1\tCAG\t3\t3\t200\t10\t1\t33.3\t3\n"
    )
    out = {k: tmp_path / f"{k}.tsv" for k in ("rare", "ex")}
    summ = tmp_path / "summary.json"
    main([
        "--p1-alleles", p1, "--p2-alleles", p2, "--plvi", str(plvi), "--min-alleles", "3", "--grid-min", "2",
        "--out-rarefaction", str(out["rare"]), "--out-examples", str(out["ex"]), "--out-summary", str(summ),
    ])
    s = json.loads(summ.read_text())
    assert s["loci"]["rarefaction_set"] == 3, "A, B, C; LOW fails the Phase 1 call rate"
    # Phase 2 alleles at A, B (REF only) and C: 6 + 1 + 2; absent from Phase 1: -3, 9, 12 at A and 30 at C.
    assert s["phase2_alleles_absent_from_phase1"] == {"alleles": 4, "of": 9, "frac": pytest.approx(4 / 9)}
    rare = list(csv.DictReader(out["rare"].open(), delimiter="\t"))
    p2_rows = [r for r in rare if r["phase"] == "phase2"]
    assert any(r["n_hap"] == "20" for r in p2_rows), "Phase 2 is evaluated at Phase 1's haplotype count"
    assert float(p2_rows[-1]["mean"]) == pytest.approx((6 + 1 + 2) / 3, rel=0.05)
    ex = list(csv.DictReader(out["ex"].open(), delimiter="\t"))
    assert {r["locus_id"] for r in ex} == {"A"}, "C is not polymorphic in Phase 1; one per motif group"
    new = [int(r["dosage_bp"]) for r in ex if r["phase"] == "phase2" and r["in_phase1"] == "0"]
    assert new == [-3, 9, 12]
