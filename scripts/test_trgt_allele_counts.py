#!/usr/bin/env python3
"""Tests for scripts/trgt_allele_counts.py (needs bcftools on PATH)."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aggregate_repeat_loci import ALLELE_COLUMNS  # noqa: E402
from trgt_allele_counts import ALLELE_COLUMNS as COUNT_COLUMNS  # noqa: E402
from trgt_allele_counts import haplotype_dosages, main, vcf_paths_from_lines  # noqa: E402

needs_bcftools = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=TRID,Number=1,Type=String,Description="TR id">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Allele length">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def _vcf(path: Path, sample: str, rows: list[tuple]) -> str:
    """rows: (trid, motifs, ref, alt, gt, al) in whatever order the VCF uses."""
    lines = [HEADER.format(sample=sample).rstrip("\n")]
    for i, (trid, motifs, ref, alt, gt, al) in enumerate(rows):
        pos = 1000 * (i + 1)
        end = pos + len(ref) - 1
        lines.append("\t".join([
            "chr1", str(pos), ".", ref, alt if alt else ".", ".", "PASS",
            f"TRID={trid};END={end};MOTIFS={motifs}", "GT:AL", f"{gt}:{al}",
        ]))
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    return str(path)


def _alleles(path: Path) -> dict[str, dict[int, int]]:
    out: dict[str, dict[int, int]] = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(r["locus_id"], {})[int(r["dosage_bp"])] = int(r["n_hap_allele"])
            out[r["locus_id"]]["_n_hap"] = int(r["n_hap"])
            out[r["locus_id"]]["_n_miss"] = int(r["n_hap_missing"])
    return out


def test_vcf_list_takes_the_path_column():
    lines = [
        "1000151\tgs://bucket/1000151.HA_v7.TR_Explorer_1.0.1.trgt.vcf.gz",
        "gs://bucket/only.vcf.gz",
        "",
        "sample_id\tgs://bucket/header-should-not-be-read.vcf.gz",
    ]
    assert vcf_paths_from_lines(lines) == [
        "gs://bucket/1000151.HA_v7.TR_Explorer_1.0.1.trgt.vcf.gz",
        "gs://bucket/only.vcf.gz",
    ]


def test_haplotype_dosages_uses_allele_length_not_gt_index():
    assert haplotype_dosages("30,36", "0/1", 30) == [0, 6]
    assert haplotype_dosages("36,36", "1|1", 30) == [6, 6]
    assert haplotype_dosages(".", "./.", 30) == [None, None]
    assert haplotype_dosages("33", "1", 30) == [3]
    assert haplotype_dosages("30,.", "0/.", 30) == [0, None]


@needs_bcftools
def test_counts_merge_and_missing_locus(tmp_path: Path):
    ref_a = "CAG" * 10
    ref_b = "AT" * 8
    s1 = _vcf(tmp_path / "s1.vcf.gz", "S1", [
        ("A", "CAG", ref_a, "CAG" * 11, "0/1", f"{len(ref_a)},{len(ref_a) + 3}"),
        ("B", "AT", ref_b, ".", "0/0", f"{len(ref_b)},{len(ref_b)}"),
    ])
    s2 = _vcf(tmp_path / "s2.vcf.gz", "S2", [
        ("A", "CAG", ref_a, "CAG" * 12, "1/1", f"{len(ref_a) + 6},{len(ref_a) + 6}"),
        ("B", "AT", ref_b, ".", "./.", "."),
    ])
    s3 = _vcf(tmp_path / "s3.vcf.gz", "S3", [
        ("A", "CAG", ref_a, ".", "./.", "."),
        ("B", "AT", ref_b, "AT" * 9, "0/1", f"{len(ref_b)},{len(ref_b) + 2}"),
    ])
    # C is in the catalog and in no VCF: two missing haplotypes per sample, not reference.
    cat = tmp_path / "catalog.bed"
    cat.write_text("chr1\t999\t1029\tA\tCAG\nchr1\t1999\t2014\tB\tAT\nchr1\t2999\t3028\tC\tCAG\n")
    sh1, sh2 = tmp_path / "a.npz", tmp_path / "b.npz"
    main(["counts", "--catalog-bed", str(cat), "--vcf", s1, "--vcf", s2, "--threads", "2", "--out", str(sh1)])
    main(["counts", "--catalog-bed", str(cat), "--vcf", s3, "--out", str(sh2)])
    out = tmp_path / "alleles.tsv"
    summary = tmp_path / "summary.json"
    main(["merge", "--catalog-bed", str(cat), "--counts", str(sh1), str(sh2), "--out", str(out), "--out-summary", str(summary)])

    header = out.read_text().splitlines()[0].split("\t")
    assert header == list(ALLELE_COLUMNS) == list(COUNT_COLUMNS)
    alleles = _alleles(out)
    assert alleles["A"] == {0: 1, 3: 1, 6: 2, "_n_hap": 4, "_n_miss": 2}
    assert alleles["B"] == {0: 3, 2: 1, "_n_hap": 4, "_n_miss": 2}
    assert alleles["C"] == {0: 0, "_n_hap": 0, "_n_miss": 6}
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert next(r for r in rows if r["locus_id"] == "A" and r["dosage_bp"] == "3")["dosage_units"] == "1"
    assert next(r for r in rows if r["locus_id"] == "B" and r["dosage_bp"] == "2")["dosage_units"] == "1"
    s = json.loads(summary.read_text())
    assert s["n_vcf"] == 3
    assert s["loci_polymorphic"] == 2
    assert s["hap_called"] == 8
    assert s["hap_missing"] == 10


@needs_bcftools
def test_off_order_and_hemizygous(tmp_path: Path):
    ref_a = "CAG" * 10
    ref_b = "AT" * 8
    # Catalog order is A then B. This VCF is B then A, and B is hemizygous.
    swapped = _vcf(tmp_path / "swapped.vcf.gz", "S2", [
        ("B", "AT", ref_b, "AT" * 9, "1", str(len(ref_b) + 2)),
        ("A", "CAG", ref_a, "CAGC" + "CAG" * 9, "1/1", f"{len(ref_a) + 1},{len(ref_a) + 1}"),
    ])
    cat = tmp_path / "catalog.bed"
    cat.write_text("chr1\t999\t1029\tA\tCAG\nchr1\t1999\t2014\tB\tAT\n")
    shard = tmp_path / "shard.npz"
    main(["counts", "--catalog-bed", str(cat), "--vcf", swapped, "--out", str(shard)])
    out = tmp_path / "alleles.tsv"
    main(["merge", "--catalog-bed", str(cat), "--counts", str(shard), "--out", str(out),
          "--out-summary", str(tmp_path / "summary.json")])
    alleles = _alleles(out)
    assert alleles["A"][1] == 2 and alleles["A"]["_n_hap"] == 2 and alleles["A"]["_n_miss"] == 0
    assert alleles["B"][2] == 1 and alleles["B"]["_n_hap"] == 1 and alleles["B"]["_n_miss"] == 0
    units = next(r["dosage_units"] for r in csv.DictReader(out.open(), delimiter="\t")
                 if r["locus_id"] == "A" and r["dosage_bp"] == "1")
    assert units == "0.3333"
