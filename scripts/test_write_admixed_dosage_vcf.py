#!/usr/bin/env python3
"""Tests for scripts/write_admixed_dosage_vcf.py.

Uses the extract-tracts-flare RU fixture; DS{k} × RU_SCALE must equal
`extract-tracts-flare --ru-baseline ref` dosages (tests/ru_dosage.rs).
"""

from __future__ import annotations

import io
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from write_admixed_dosage_vcf import write_admixed_vcf  # noqa: E402

FIXTURE = ROOT / "tractor_mix" / "extract_tracts_flare" / "testdata" / "ru_dosage.vcf"


def _parse(text: str):
    header, rows = None, {}
    for line in text.splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#CHROM"):
            header = f
            continue
        scale = float(f[7].split("RU_SCALE=")[1])
        keys = f[8].split(":")
        per_sample = []
        for cell in f[9:]:
            vals = dict(zip(keys, cell.split(":")))
            per_sample.append(vals)
        rows[f[2]] = (scale, keys, per_sample)
    return header, rows


def _unscaled(row, k):
    scale, _, per_sample = row
    return [round(float(s[f"DS{k}"]) * scale, 5) for s in per_sample]


def test_matches_extract_ref_baseline():
    buf = io.StringIO()
    stats = write_admixed_vcf(str(FIXTURE), buf, num_ancs=3)
    header, rows = _parse(buf.getvalue())
    assert header[9:] == ["S1", "S2", "S3"]
    assert "snv1" not in rows, "non-RU_TEST records are skipped"
    assert stats == {"records_in": 5, "ru_records": 4, "written": 4, "monomorphic_skipped": 0}
    assert rows["vntr_ins"][1] == ["DS1", "DS2", "DS3", "ANC1", "ANC2", "ANC3"]
    # Same numbers as tests/ru_dosage.rs ru_dosage_ref_baseline_is_default.
    assert _unscaled(rows["vntr_ins"], 1) == [0, 6, 0]
    assert _unscaled(rows["vntr_ins"], 2) == [3, 0, 0]
    assert _unscaled(rows["vntr_del"], 1) == [0, 0, -5]
    assert _unscaled(rows["vntr_del"], 2) == [-5, 0, -5]
    assert _unscaled(rows["vntr_del"], 3) == [0, 0, 0]
    assert _unscaled(rows["vntr_multi"], 1) == [1, 0, 0]
    assert _unscaled(rows["vntr_multi"], 2) == [3, 1, 0]
    assert _unscaled(rows["rel_ins"], 1) == [1, 0, 0]
    assert _unscaled(rows["rel_ins"], 2) == [0, 2, 0]
    # Scale = max |x| at the locus; scaled |DS| never exceeds ANC.
    assert rows["vntr_ins"][0] == 3 and rows["vntr_del"][0] == 5 and rows["vntr_multi"][0] == 3
    for scale, _, per_sample in rows.values():
        for s in per_sample:
            for k in (1, 2, 3):
                assert abs(float(s[f"DS{k}"])) <= float(s[f"ANC{k}"]) + 1e-12
    # Haplotype counts per ancestry.
    s1, s2, s3 = rows["vntr_ins"][2]
    assert (s1["ANC1"], s1["ANC2"], s1["ANC3"]) == ("1", "1", "0")
    assert (s2["ANC1"], s3["ANC2"]) == ("2", "2")


def test_keep_list_reorders_and_errors():
    buf = io.StringIO()
    write_admixed_vcf(str(FIXTURE), buf, num_ancs=3, keep=["S3", "S1"])
    header, rows = _parse(buf.getvalue())
    assert header[9:] == ["S3", "S1"]
    assert _unscaled(rows["vntr_ins"], 1) == [0, 0]
    try:
        write_admixed_vcf(str(FIXTURE), io.StringIO(), num_ancs=3, keep=["nope"])
    except SystemExit as exc:
        assert "not in VCF" in str(exc)
    else:
        raise AssertionError("expected missing-sample error")


def test_rejects_unphased_and_bad_ancestry():
    text = FIXTURE.read_text()
    with tempfile.TemporaryDirectory() as td:
        bad = Path(td) / "bad.vcf"
        bad.write_text(text.replace("0|1:0:1\t1|1:0:0\t0|0:1:1", "0/1:0:1\t1|1:0:0\t0|0:1:1"))
        try:
            write_admixed_vcf(str(bad), io.StringIO(), num_ancs=3)
        except ValueError as exc:
            assert "unphased" in str(exc)
        else:
            raise AssertionError("expected unphased error")
        try:
            write_admixed_vcf(str(FIXTURE), io.StringIO(), num_ancs=2)
        except ValueError as exc:
            assert "ancestry 2 outside" in str(exc)
        else:
            raise AssertionError("expected ancestry range error")


if __name__ == "__main__":
    test_matches_extract_ref_baseline()
    test_keep_list_reorders_and_errors()
    test_rejects_unphased_and_bad_ancestry()
    print("ok")
