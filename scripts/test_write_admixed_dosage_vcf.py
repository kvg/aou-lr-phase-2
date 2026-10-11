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
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        scale = float(info["RU_SCALE"])
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
    assert stats == {
        "records_in": 5, "ru_records": 4, "written": 4, "monomorphic_skipped": 0,
        "missingness_skipped": 0, "hap_missing": 0, "hap_total": 24, "coding": "ru",
    }
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
            write_admixed_vcf(str(bad), io.StringIO(), num_ancs=3, missing="error")
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


def _missing_fixture(td: Path) -> Path:
    """One locus; S1 half-missing, S2 unphased-ambiguous, S3 fully called."""
    vcf = td / "missing.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr22,length=100000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "chr22\t900\tlocus1\tN\t<TR:+6bp>,<TR:-3bp>\t.\tPASS\tRU_TEST;RU_DOSAGE=2,-1\t"
        "GT:AN1:AN2\t1|.:0:1\t./.:0:0\t1|2:1:1\t0|0:0:0\n"
    )
    return vcf


def test_missing_ref_zeroes_dosage_and_shrinks_anc():
    with tempfile.TemporaryDirectory() as td:
        vcf = _missing_fixture(Path(td))
        buf = io.StringIO()
        stats = write_admixed_vcf(str(vcf), buf, num_ancs=2, max_hap_missing=0.5)
    _, rows = _parse(buf.getvalue())
    scale, _, per_sample = rows["locus1"]
    s1, s2, s3, s4 = per_sample
    # S1: hap1 = +2 units on ancestry 0; hap2 missing, so it leaves ANC2 as well.
    assert _unscaled(rows["locus1"], 1) == [2, 0, 0, 0]
    assert (s1["ANC1"], s1["ANC2"]) == ("1", "0")
    # S2 is unphased-ambiguous: no haplotype enters either ancestry.
    assert (s2["DS1"], s2["DS2"], s2["ANC1"], s2["ANC2"]) == ("0", "0", "0", "0")
    # S3 is fully called on ancestry 1: +2 and -1 units sum to +1.
    assert _unscaled(rows["locus1"], 2) == [0, 0, 1, 0]
    assert (s3["ANC1"], s3["ANC2"]) == ("0", "2")
    # S4 is reference on both haplotypes and still contributes its denominator.
    assert (s4["ANC1"], s4["ANC2"]) == ("2", "0")
    assert stats["hap_missing"] == 3 and stats["hap_total"] == 8
    assert stats["written"] == 1 and stats["missingness_skipped"] == 0


def test_missingness_cap_skips_locus_and_is_recorded():
    with tempfile.TemporaryDirectory() as td:
        vcf = _missing_fixture(Path(td))
        buf = io.StringIO()
        stats = write_admixed_vcf(str(vcf), buf, num_ancs=2, max_hap_missing=0.2)
        assert stats["missingness_skipped"] == 1 and stats["written"] == 0
        assert "locus1" not in buf.getvalue()

        buf = io.StringIO()
        write_admixed_vcf(str(vcf), buf, num_ancs=2, max_hap_missing=0.5)
        info = [ln for ln in buf.getvalue().splitlines() if ln.startswith("chr22")][0].split("\t")[7]
        assert "RU_HAP_MISS=0.375" in info, info

        try:
            write_admixed_vcf(str(vcf), io.StringIO(), num_ancs=2, missing="error")
        except ValueError as exc:
            assert "unphased" in str(exc) or "missing allele" in str(exc)
        else:
            raise AssertionError("expected strict-mode refusal")


def test_sign_per_alt_and_signed_locus_dosage():
    head, body = FIXTURE.read_text().split("#CHROM", 1)
    cols = "#CHROM" + body.splitlines()[0]
    records = [
        # Sequence-resolved deletion with no SVTYPE: shorter ALT is a loss.
        "chr22\t800\tsmall_del\tACAGCAG\tA\t.\tPASS\tRU_TEST;PERIOD=3;RU=2\tGT:AN1:AN2\t0|1:0:1\t0|0:0:0\t1|1:1:1",
        # Aggregated locus: signed per-ALT dosage, one gain and one loss.
        "chr22\t900\tlocus1\tN\t<TR:+6>,<TR:-3>\t.\tPASS\tRU_TEST;RU_DOSAGE=2,-1\tGT:AN1:AN2\t1|2:0:1\t0|2:1:1\t0|0:0:0",
    ]
    with tempfile.TemporaryDirectory() as td:
        vcf = Path(td) / "signed.vcf"
        vcf.write_text(head + cols + "\n" + "\n".join(records) + "\n")
        buf = io.StringIO()
        write_admixed_vcf(str(vcf), buf, num_ancs=3)
    _, rows = _parse(buf.getvalue())
    assert _unscaled(rows["small_del"], 1) == [0, 0, 0]
    assert _unscaled(rows["small_del"], 2) == [-2, 0, -4]
    assert _unscaled(rows["locus1"], 1) == [2, 0, 0]
    assert _unscaled(rows["locus1"], 2) == [-1, -1, 0]
    assert rows["locus1"][0] == 2


def test_carrier_coding_counts_non_reference_haplotypes():
    buf = io.StringIO()
    stats = write_admixed_vcf(str(FIXTURE), buf, num_ancs=3, coding="carrier")
    assert stats["coding"] == "carrier"
    _, rows = _parse(buf.getvalue())
    # RU_SCALE is 1, so DS{k} is the raw non-reference haplotype count.
    for scale, _, per_sample in rows.values():
        assert scale == 1.0
        for s in per_sample:
            for k in (1, 2, 3):
                ds, anc = float(s[f"DS{k}"]), float(s[f"ANC{k}"])
                assert ds == int(ds) and 0 <= ds <= anc <= 2
    # vntr_ins under ru coding: S2 carries 6 units on anc1, S1 3 units on anc2.
    assert _unscaled(rows["vntr_ins"], 1) == [0, 2, 0]
    assert _unscaled(rows["vntr_ins"], 2) == [1, 0, 0]
    # Deletions are carriers too: the sign is dropped, not the haplotype.
    assert _unscaled(rows["vntr_del"], 1) == [0, 0, 1]
    assert _unscaled(rows["vntr_del"], 2) == [1, 0, 1]
    # This is the information the dosage arm adds: on ancestry 2, S1 carries a
    # 3-unit allele and S2 a 1-unit allele (ru: 1 vs 1/3), but both are one
    # carrier haplotype, so the biallelic arm cannot tell them apart.
    assert _unscaled(rows["vntr_multi"], 2) == [1, 1, 0]
    info = [ln for ln in buf.getvalue().splitlines() if not ln.startswith("#")][0].split("\t")[7]
    assert "RU_CODING=carrier" in info, info


def test_codings_share_loci_and_haplotypes():
    """The 2x2 needs the arms matched: same loci, same ANC{k}, same carriers.

    One-directional on DS: a non-zero repeat dosage implies a carrier, but not
    the reverse — see test_carrier_coding_differs_only_by_cancellation.
    """
    out = {}
    for coding in ("ru", "carrier"):
        buf = io.StringIO()
        write_admixed_vcf(str(FIXTURE), buf, num_ancs=3, coding=coding)
        out[coding] = _parse(buf.getvalue())[1]
    assert list(out["ru"]) == list(out["carrier"])
    for locus in out["ru"]:
        ru_rows = out["ru"][locus][2]
        ca_rows = out["carrier"][locus][2]
        for ru_cell, ca_cell in zip(ru_rows, ca_rows):
            for k in (1, 2, 3):
                assert ru_cell[f"ANC{k}"] == ca_cell[f"ANC{k}"]
                if float(ru_cell[f"DS{k}"]) != 0:
                    assert float(ca_cell[f"DS{k}"]) != 0


def test_carrier_coding_differs_only_by_cancellation():
    """A gain and a loss of equal size on one ancestry cancel under ru coding.

    FELIX's carrier QC counts samples with a non-zero DS{k}, so this sample
    enters the biallelic arm's QC denominator and not the dosage arm's. It is
    the one way the two arms' carrier counts can disagree, and the reason
    compare_ru_arms.py reports n_loci_*_only instead of assuming a match.
    """
    with tempfile.TemporaryDirectory() as td:
        vcf = Path(td) / "cancel.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr22,length=100000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr22\t900\tlocus1\tN\t<TR:+2>,<TR:-2>\t.\tPASS\tRU_TEST;RU_DOSAGE=2,-2\t"
            "GT:AN1:AN2\t1|2:0:0\t1|0:0:0\n"
        )
        got = {}
        for coding in ("ru", "carrier"):
            buf = io.StringIO()
            write_admixed_vcf(str(vcf), buf, num_ancs=2, coding=coding)
            got[coding] = _parse(buf.getvalue())[1]["locus1"][2]
    (ru_s1, ru_s2), (ca_s1, ca_s2) = got["ru"], got["carrier"]
    # S1 carries both alleles on ancestry 1: net zero units, two carrier haplotypes.
    assert ru_s1["DS1"] == "0" and ca_s1["DS1"] == "2"
    assert ru_s1["ANC1"] == ca_s1["ANC1"] == "2"
    # S2 is unambiguous in both arms.
    assert float(ru_s2["DS1"]) != 0 and float(ca_s2["DS1"]) != 0


def test_carrier_coding_keeps_missingness_policy():
    with tempfile.TemporaryDirectory() as td:
        vcf = _missing_fixture(Path(td))
        buf = io.StringIO()
        stats = write_admixed_vcf(str(vcf), buf, num_ancs=2, max_hap_missing=0.5, coding="carrier")
    _, rows = _parse(buf.getvalue())
    s1, s2, s3, s4 = rows["locus1"][2]
    assert (s1["DS1"], s1["ANC1"], s1["ANC2"]) == ("1", "1", "0")
    assert (s2["DS1"], s2["DS2"], s2["ANC1"], s2["ANC2"]) == ("0", "0", "0", "0")
    # S3 carries +2 and -1 units on ancestry 2: one signed unit, but two carriers.
    assert (s3["DS2"], s3["ANC2"]) == ("2", "2")
    assert (s4["DS1"], s4["ANC1"]) == ("0", "2")
    assert stats["hap_missing"] == 3 and stats["written"] == 1


def test_rejects_unknown_coding():
    try:
        write_admixed_vcf(str(FIXTURE), io.StringIO(), num_ancs=3, coding="biallelic")
    except ValueError as exc:
        assert "coding must be" in str(exc)
    else:
        raise AssertionError("expected coding validation error")


if __name__ == "__main__":
    test_matches_extract_ref_baseline()
    test_keep_list_reorders_and_errors()
    test_rejects_unphased_and_bad_ancestry()
    test_missing_ref_zeroes_dosage_and_shrinks_anc()
    test_missingness_cap_skips_locus_and_is_recorded()
    test_sign_per_alt_and_signed_locus_dosage()
    test_carrier_coding_counts_non_reference_haplotypes()
    test_codings_share_loci_and_haplotypes()
    test_carrier_coding_differs_only_by_cancellation()
    test_carrier_coding_keeps_missingness_policy()
    test_rejects_unknown_coding()
    print("ok")
