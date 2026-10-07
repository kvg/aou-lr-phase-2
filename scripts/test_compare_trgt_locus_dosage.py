#!/usr/bin/env python3
"""Tests for scripts/compare_trgt_locus_dosage.py (needs bcftools on PATH)."""

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

from aggregate_repeat_loci import main as aggregate  # noqa: E402
from compare_trgt_locus_dosage import main as compare  # noqa: E402
from test_aggregate_repeat_loci import _write_vcf  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

TRGT_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
##INFO=<ID=TRID,Number=1,Type=String,Description="TR id">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def _trgt(path: Path, sample: str, calls: dict[str, tuple[int, int]]) -> str:
    starts = {"L1": 500, "L2": 3001, "L3": 8001}
    lines = [TRGT_HEADER.format(sample=sample).rstrip("\n")]
    for lid, (d1, d2) in calls.items():
        ref = "C" * 30
        alts, gt = [], []
        for d in (d1, d2):
            if d == 0:
                gt.append("0")
                continue
            seq = "C" * (30 + d)
            if seq not in alts:
                alts.append(seq)
            gt.append(str(alts.index(seq) + 1))
        lines.append("\t".join([
            "chr22", str(starts[lid]), ".", ref, ",".join(alts) or ".", ".", "PASS",
            f"TRID={lid};END={starts[lid] + 29};MOTIFS=CAG", "GT", "/".join(gt),
        ]))
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(path)], check=True)
    return str(path)


def test_concordance(tmp_path: Path):
    vcf = _write_vcf(tmp_path / "in.vcf.gz")
    bed = tmp_path / "catalog.bed"
    bed.write_text("chr22\t499\t530\tL1\tCAG\nchr22\t3000\t3020\tL2\tAT\nchr22\t8000\t8030\tL3\tCAG\n")
    loci_vcf = tmp_path / "loci.vcf"
    aggregate([
        "--vcf", str(vcf), "--catalog-bed", str(bed), "--chrom", "chr22", "--exclude-prefixes", "HG",
        "--out-alleles", str(tmp_path / "a.tsv"), "--out-vcf", str(loci_vcf), "--out-summary", str(tmp_path / "s.json"),
    ])
    # Integrated: L1 S1 63|60 S2 3|-6 S3 3|0; L2 S1 ambiguous S2 0|2 S3 -4|-4; L3 not in the locus VCF.
    trgt = {
        "S1": {"L1": (3, 60), "L2": (0, 0)},
        "S2": {"L1": (-6, 3), "L2": (2, 0), "L3": (0, 0)},
        "S3": {"L1": (0, 3), "L2": (-4, -2)},
    }
    table = tmp_path / "trgt.tsv"
    table.write_text("".join(f"{s}\t{_trgt(tmp_path / f'{s}.trgt.vcf.gz', s, c)}\n" for s, c in trgt.items()))
    out_loci, out_sum = tmp_path / "loci.tsv", tmp_path / "summary.json"
    compare([
        "--locus-vcf", str(loci_vcf), "--catalog-bed", str(bed), "--trgt-tsv", str(table), "--chrom", "chr22",
        "--out-loci", str(out_loci), "--out-summary", str(out_sum),
    ])
    rows = {r["locus_id"]: r for r in csv.DictReader(out_loci.open(), delimiter="\t")}
    assert (rows["L1"]["n_pairs"], rows["L1"]["n_exact"], rows["L1"]["n_overcount"]) == ("3", "2", "1")
    assert (rows["L2"]["n_pairs"], rows["L2"]["n_exact"], rows["L2"]["n_within_unit"]) == ("2", "1", "2")
    assert rows["L3"]["n_records"] == "0" and rows["L3"]["n_exact"] == "1"
    s = json.loads(out_sum.read_text())
    assert s["n_samples"] == 3 and s["overall"]["pairs"] == 6
    assert s["by_integrated_records"]["3+"]["overcount"] == pytest.approx(1 / 3)
