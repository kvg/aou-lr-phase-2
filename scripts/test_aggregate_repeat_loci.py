#!/usr/bin/env python3
"""Tests for scripts/aggregate_repeat_loci.py (needs bcftools on PATH)."""

from __future__ import annotations

import csv
import gzip
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aggregate_repeat_loci import allele_bp, main  # noqa: E402
from write_admixed_dosage_vcf import write_admixed_vcf  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=RU_TEST,Number=0,Type=Flag,Description="RU">
##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Period">
##INFO=<ID=MOTIF,Number=1,Type=String,Description="Motif">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry hap 1">
##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry hap 2">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tHG00001
"""

# Locus L1 = chr22:[499, 530) CAG; locus L2 = chr22:[3000, 3020) AT.
# Haplotype sums (bp): L1  S1 63|60  S2 3|-6  S3 3|0  S4 missing (r3 is ./.)
#                      L2  S1 ambiguous (unphased het at two records)  S2 0|2  S3 -4|-4  S4 0|0
RECORDS = [
    ("500", "r1", "A", "ACAG", ".", ["1|0:0:1", "0|0:1:1", "1|1:2:2", "0|1:0:0", "0|0:0:0"]),
    ("505", "r2", "ACAGCAG", "A", ".", ["0|0:0:1", "0|1:1:1", "0|0:2:2", "0|0:0:0", "0|0:0:0"]),
    ("510", "r3", "N", "<INS>", "SVTYPE=INS;SVLEN=60", ["1|1:0:1", "0|0:1:1", "0|0:2:2", "./.:0:0", "0|0:0:0"]),
    ("512", "snv", "C", "T", ".", ["0|1:0:1", "0|0:1:1", "0|0:2:2", "0|0:0:0", "1|1:0:0"]),
    ("515", "r5", "ACAG", "A,ACAGCAG", ".", ["0|0:0:1", "2|0:1:1", "0|1:2:2", "0|0:0:0", "0|0:0:0"]),
    ("2000", "lone", "A", "ATT", "RU_TEST;PERIOD=1;MOTIF=T", ["0|1:0:1", "0|0:1:1", "0|0:2:2", "0|0:0:0", "0|0:0:0"]),
    ("3005", "r6", "TATAT", "T", ".", ["0/1:0:1", "0|0:1:1", "1|1:2:2", "0|0:0:0", "0|0:0:0"]),
    ("3010", "r7", "T", "TAT", ".", ["0/1:0:1", "0|1:1:1", "0|0:2:2", "0|0:0:0", "0|0:0:0"]),
]


def _write_vcf(path: Path) -> Path:
    lines = [HEADER.rstrip("\n")]
    for pos, rid, ref, alt, info, gts in RECORDS:
        lines.append("\t".join(["chr22", pos, rid, ref, alt, ".", "PASS", info, "GT:AN1:AN2", *gts]))
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(path)], check=True)
    return path


def _run(tmp_path: Path, *extra: str):
    vcf = _write_vcf(tmp_path / "in.vcf.gz")
    bed = tmp_path / "catalog.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\nchr22\t3000\t3020\tL2\tAT\n")
    alleles = tmp_path / "alleles.tsv.gz"
    out_vcf = tmp_path / "loci.vcf"
    summary = tmp_path / "summary.json"
    main([
        "--vcf", str(vcf), "--catalog-bed", str(bed), "--chrom", "chr22", "--exclude-prefixes", "HG",
        "--out-alleles", str(alleles), "--out-vcf", str(out_vcf), "--out-summary", str(summary), *extra,
    ])
    with gzip.open(alleles, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows, out_vcf.read_text(), json.loads(summary.read_text())


def test_symbolic_deletion_span_uses_svlen_when_end_is_absent():
    from aggregate_repeat_loci import record_span
    bp = allele_bp("N", ["<DEL>"], "DEL", "-400", ".", 1000)
    assert bp.tolist() == [0, -400]
    assert record_span(1000, "N", bp, ".") == (999, 1399)


def test_allele_bp_signs():
    assert allele_bp("ACAG", ["A", "ACAGCAG"], ".", ".", ".", 1).tolist() == [0, -3, 3]
    assert allele_bp("N", ["<DEL>"], "DEL", "-10", "1010", 1000).tolist() == [0, -10]
    assert allele_bp("N", ["<DEL>"], "DEL", ".", "1010", 1000).tolist() == [0, -10]
    assert allele_bp("N", ["<DUP>"], ".", "40", ".", 1).tolist() == [0, 40]
    assert allele_bp("N", ["<INV>"], "INV", "500", ".", 1) is None
    assert allele_bp("C", ["T"], ".", ".", ".", 1) is None
    assert allele_bp("ACAG", ["A", "*"], ".", ".", ".", 1).tolist() == [0, -3, 0]


def test_locus_sums_alleles_and_diagnostics(tmp_path: Path):
    rows, _vcf, summary = _run(tmp_path)
    by_locus: dict[str, dict[int, int]] = {}
    for r in rows:
        by_locus.setdefault(r["locus_id"], {})[int(r["dosage_bp"])] = int(r["n_hap_allele"])
    assert by_locus["L1"] == {-6: 1, 0: 1, 3: 2, 60: 1, 63: 1}
    assert by_locus["L2"] == {-4: 2, 0: 3, 2: 1}
    l1 = next(r for r in rows if r["locus_id"] == "L1")
    assert (l1["n_hap"], l1["n_hap_missing"], l1["n_records"], l1["period"]) == ("6", "2", "4", "3")
    assert next(r for r in rows if r["locus_id"] == "L1" and r["dosage_bp"] == "63")["dosage_units"] == "21"
    c = summary["counts"]
    assert summary["n_samples"] == 4, "HG00001 is excluded"
    # The SNV is dropped by the bcftools prefilter and the record outside every locus is never fetched.
    assert c["records_read"] == 6 and c["records_assigned"] == 6 and c["records_unassigned"] == 0
    assert c["loci_with_records"] == 2 and c["loci_polymorphic"] == 2
    assert c["hap_unphased_ambiguous"] == 2
    assert c["hap_multi_record"] == 2
    assert c["hap_possible_duplicate"] == 1


def test_unphased_hom_plus_het_resolves(tmp_path: Path, monkeypatch):
    # S1: r6 1/1 (-4 on both) + r7 0/1 (+2) -> -2|-4, no ambiguity; S2: r6 0/1 + r7 0/1 -> ambiguous.
    import test_aggregate_repeat_loci as m
    records = [
        ("3005", "r6", "TATAT", "T", ".", ["1/1:0:1", "0/1:1:1", "0/0:2:2", "0/0:0:0", "0/0:0:0"]),
        ("3010", "r7", "T", "TAT", ".", ["0/1:0:1", "0/1:1:1", "0/0:2:2", "0/0:0:0", "0/0:0:0"]),
    ]
    monkeypatch.setattr(m, "RECORDS", records)
    rows, _vcf, summary = _run(tmp_path)
    l2 = {int(r["dosage_bp"]): int(r["n_hap_allele"]) for r in rows if r["locus_id"] == "L2"}
    assert l2 == {-4: 1, -2: 1, 0: 4}
    assert summary["counts"]["hap_unphased_ambiguous"] == 2


def test_ignore_phase(tmp_path: Path):
    # S3 is phased 0|1 at r5 (L1) with a 1|1 at r1: still resolvable. S1 at L1: r1 1|0 + r3 1|1 -> resolvable.
    # S2 at L1: r2 0|1 + r5 2|0, two hets -> ambiguous once phase is ignored.
    _rows, _vcf, summary = _run(tmp_path, "--ignore-phase")
    assert summary["ignore_phase"] is True
    assert summary["counts"]["hap_unphased_ambiguous"] == 2 + 2


def _alleles(tmp_path: Path, name: str, *vcf_args: str) -> tuple[dict, dict]:
    bed = tmp_path / "catalog.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\nchr22\t3000\t3020\tL2\tAT\n")
    out, summ = tmp_path / f"{name}.tsv", tmp_path / f"{name}.json"
    main([*vcf_args, "--catalog-bed", str(bed), "--chrom", "chr22", "--exclude-prefixes", "HG",
          "--apply-filters", "PASS,.", "--out-alleles", str(out), "--out-summary", str(summ)])
    got: dict[str, dict[int, int]] = {}
    for r in csv.DictReader(out.open(), delimiter="\t"):
        got.setdefault(r["locus_id"], {})[int(r["dosage_bp"])] = int(r["n_hap_allele"])
    return got, json.loads(summ.read_text())


def test_small_and_sv_sources_count_each_allele_once(tmp_path: Path):
    vcf = str(_write_vcf(tmp_path / "in.vcf.gz"))
    whole, _ = _alleles(tmp_path, "whole", "--vcf", vcf)
    both, s_both = _alleles(tmp_path, "both", "--small-vcf", vcf, "--sv-vcf", vcf)
    assert both == whole, "each allele is kept from exactly one of the two copies"
    assert s_both["counts"]["alleles_out_of_band"] == 6 + 1, "SV copy drops 6 small ALTs, small copy drops the 60 bp"
    small, _ = _alleles(tmp_path, "small", "--small-vcf", vcf)
    assert 60 not in small["L1"] and 63 not in small["L1"]
    sv, _ = _alleles(tmp_path, "sv", "--sv-vcf", vcf)
    assert sv["L1"] == {0: 4, 60: 2}, "only S1's 60 bp insertion; S4 is missing at r3"


def test_gvcf_reference_block_is_not_read(tmp_path: Path):
    header = HEADER.replace("S1\tS2\tS3\tS4\tHG00001", "S1")
    lines = [header.rstrip("\n")]
    ref = "N" * 5000
    lines.append("\t".join(["chr22", "500", "block", ref, "<NON_REF>", ".", "PASS", "END=5500", "GT", "0/0"]))
    lines.append("\t".join(["chr22", "510", "ins", "A", "ACAG", ".", "PASS", ".", "GT", "0/1"]))
    plain = tmp_path / "in.vcf"
    plain.write_text("\n".join(lines) + "\n")
    vcf = tmp_path / "in.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(vcf)], check=True)
    bed = tmp_path / "catalog.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\n")
    summary = tmp_path / "summary.json"
    main([
        "--vcf", str(vcf), "--catalog-bed", str(bed), "--chrom", "chr22",
        "--out-alleles", str(tmp_path / "alleles.tsv"), "--out-summary", str(summary),
    ])
    counts = json.loads(summary.read_text())["counts"]
    assert counts["records_read"] == 1
    assert counts["records_length_change"] == 1


def test_gvcf_without_info_end_still_reads(tmp_path: Path):
    header = "\n".join(
        ln for ln in HEADER.replace("S1\tS2\tS3\tS4\tHG00001", "S1").splitlines()
        if not ln.startswith("##INFO=<ID=SVTYPE") and not ln.startswith("##INFO=<ID=SVLEN") and not ln.startswith("##INFO=<ID=END")
    )
    plain = tmp_path / "in.vcf"
    plain.write_text(header + "\n" + "\t".join(["chr22", "510", "ins", "A", "ACAG", ".", "PASS", ".", "GT", "0/1"]) + "\n")
    vcf = tmp_path / "in.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(vcf)], check=True)
    bed = tmp_path / "catalog.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\n")
    summary = tmp_path / "summary.json"
    main([
        "--vcf", str(vcf), "--catalog-bed", str(bed), "--chrom", "chr22",
        "--out-alleles", str(tmp_path / "alleles.tsv"), "--out-summary", str(summary),
    ])
    assert json.loads(summary.read_text())["counts"]["records_length_change"] == 1


def test_keep_uncatalogued_ru_record(tmp_path: Path):
    rows, vcf, summary = _run(tmp_path, "--keep-uncatalogued-ru")
    c = summary["counts"]
    assert c["records_read"] == 7 and c["records_uncatalogued_ru"] == 1 and c["loci_with_records"] == 3
    lone = {int(r["dosage_bp"]): r for r in rows if r["locus_id"] == "lone"}
    assert lone[2]["dosage_units"] == "2" and lone[2]["period"] == "1" and lone[2]["n_hap_allele"] == "1"
    ids = [ln.split("\t")[2] for ln in vcf.splitlines() if not ln.startswith("#")]
    assert ids == ["L1", "lone", "L2"], "locus VCF stays in position order"


def test_locus_vcf_feeds_dosage_writer(tmp_path: Path):
    _rows, vcf, _s = _run(tmp_path, "--with-ancestry", "--missing-as-ref")
    recs = {ln.split("\t")[2]: ln.split("\t") for ln in vcf.splitlines() if not ln.startswith("#")}
    l1 = recs["L1"]
    assert l1[4] == "<TR:-6bp>,<TR:+3bp>,<TR:+60bp>,<TR:+63bp>"
    assert "RU_DOSAGE=-2,1,20,21" in l1[7]
    assert l1[8] == "GT:AN1:AN2"
    # --missing-as-ref: S4's ./. at r3 counts as REF, so S4 is 0|+3.
    assert l1[9:] == ["4|3:0:1", "2|1:1:1", "2|0:2:2", "0|2:0:0"]
    assert recs["L2"][9] == "./.:0:1", "ambiguous unphased sample stays missing"
    phased_only = "\n".join(ln for ln in vcf.splitlines() if not ln.startswith("chr22\t3001")) + "\n"
    path = tmp_path / "phased.vcf"
    path.write_text(phased_only)
    buf = io.StringIO()
    stats = write_admixed_vcf(str(path), buf, num_ancs=3)
    assert stats["written"] == 1
    cells = [ln for ln in buf.getvalue().splitlines() if not ln.startswith("#")][0].split("\t")
    assert cells[7] == "RU_SCALE=21;RU_HAP_MISS=0"
    s1 = dict(zip(cells[8].split(":"), cells[9].split(":")))
    assert round(float(s1["DS1"]) * 21, 6) == 21 and round(float(s1["DS2"]) * 21, 6) == 20


def _two_digit_vcf(tmp_path: Path) -> Path:
    """Two records at one locus; ancestry codes 12 and 7, inconsistent at S2 hap1."""
    rows = [
        ("500", "r1", "A", "ACAG", ".", ["1|0:12:7", "0|0:7:7", "0|0:12:12"]),
        ("505", "r2", "A", "ACAGCAG", ".", ["0|1:12:7", "0|0:12:7", "0|0:12:12"]),
    ]
    header = "\n".join(
        ln for ln in HEADER.strip().splitlines()
        if not ln.startswith("#CHROM")
    ) + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
    lines = [header.rstrip("\n")]
    for pos, rid, ref, alt, info, gts in rows:
        lines.append("\t".join(["chr22", pos, rid, ref, alt, ".", "PASS", info, "GT:AN1:AN2", *gts]))
    plain = tmp_path / "wide.vcf"
    plain.write_text("\n".join(lines) + "\n")
    gz = tmp_path / "wide.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(gz), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(gz)], check=True)
    return gz


def _run_wide(tmp_path: Path, *extra: str):
    gz = _two_digit_vcf(tmp_path)
    bed = tmp_path / "wide.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\n")
    out_vcf = tmp_path / "wide.loci.vcf"
    summary = tmp_path / "wide.summary.json"
    rc = main([
        "--vcf", str(gz), "--catalog-bed", str(bed), "--chrom", "chr22",
        "--with-ancestry", "--out-vcf", str(out_vcf),
        "--out-alleles", str(tmp_path / "wide.alleles.tsv"),
        "--out-summary", str(summary), *extra,
    ])
    assert rc == 0
    return out_vcf.read_text(), json.loads(summary.read_text())


def test_two_digit_ancestry_codes_round_trip(tmp_path: Path):
    vcf, summary = _run_wide(tmp_path)
    rec = [ln.split("\t") for ln in vcf.splitlines() if ln.startswith("chr22")][0]
    assert rec[8] == "GT:AN1:AN2"
    # Codes 12 and 7 survive; reading one byte would have given 1 and 7.
    assert rec[9].endswith(":12:7"), rec[9]
    assert summary["anc_code_max"] == 12
    assert summary["counts"]["anc_code_max"] == 12
    # The writer accepts them and puts the dosage on the right ancestry.
    path = tmp_path / "wide.for_writer.vcf"
    path.write_text(vcf)
    buf = io.StringIO()
    stats = write_admixed_vcf(str(path), buf, num_ancs=13)
    assert stats["written"] == 1
    cells = [ln for ln in buf.getvalue().splitlines() if not ln.startswith("#")][0].split("\t")
    s1 = dict(zip(cells[8].split(":"), cells[9].split(":")))
    scale = float(dict(kv.split("=", 1) for kv in cells[7].split(";") if "=" in kv)["RU_SCALE"])
    # S1 carries +3 bp (1 unit) on hap1 (ancestry 12) and +6 bp (2 units) on hap2 (ancestry 7).
    assert round(float(s1["DS13"]) * scale, 6) == 1.0
    assert round(float(s1["DS8"]) * scale, 6) == 2.0
    assert (s1["ANC13"], s1["ANC8"]) == ("1", "1")


def test_ancestry_inconsistency_is_counted(tmp_path: Path):
    _vcf, summary = _run_wide(tmp_path)
    # S2 hap1 is ancestry 7 at r1 and 12 at r2; one haplotype disagrees.
    assert summary["counts"]["hap_ancestry_inconsistent"] == 1
    assert summary["hap_ancestry_inconsistent_frac"] == pytest.approx(1 / 6)


def test_num_ancs_rejects_out_of_range_code(tmp_path: Path):
    with pytest.raises(SystemExit) as exc:
        _run_wide(tmp_path, "--num-ancs", "8")
    assert "ancestry code 12" in str(exc.value)


def test_ambiguous_locus_survives_the_writer(tmp_path: Path):
    """The ./. locus the aggregator emits is tested, not refused (--missing ref)."""
    _rows, vcf, _s = _run(tmp_path, "--with-ancestry")
    path = tmp_path / "all.vcf"
    path.write_text(vcf if vcf.endswith("\n") else vcf + "\n")
    buf = io.StringIO()
    stats = write_admixed_vcf(str(path), buf, num_ancs=3, max_hap_missing=0.5)
    assert stats["written"] == 2, "both loci tested, including the ambiguous one"
    assert stats["hap_missing"] > 0 and stats["missingness_skipped"] == 0
    recs = {ln.split("\t")[2]: ln.split("\t") for ln in buf.getvalue().splitlines()
            if not ln.startswith("#")}
    l2 = recs["L2"]
    s1 = dict(zip(l2[8].split(":"), l2[9].split(":")))
    # The ambiguous sample contributes no dosage and no denominator.
    assert {s1["DS1"], s1["DS2"], s1["DS3"]} == {"0"}
    assert (s1["ANC1"], s1["ANC2"], s1["ANC3"]) == ("0", "0", "0")

    # Strict mode still refuses it, so the old behaviour stays reachable.
    try:
        write_admixed_vcf(str(path), io.StringIO(), num_ancs=3, missing="error")
    except ValueError as exc:
        assert "unphased" in str(exc) or "missing allele" in str(exc)
    else:
        raise AssertionError("expected strict-mode refusal")
