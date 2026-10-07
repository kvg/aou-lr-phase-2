#!/usr/bin/env python3
"""Tests for scripts/trgt_plvi.py (needs bcftools on PATH)."""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from trgt_plvi import lps_bp, main  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=TRID,Number=1,Type=String,Description="TR id">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Allele length">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

CAG = "CAG" * 10
INTERRUPTED = "CAG" * 4 + "CAA" + "CAG" * 6


def _vcf(path: Path, sample: str, rows: list[tuple[str, str, str, str]]) -> str:
    lines = [HEADER.format(sample=sample).rstrip("\n")]
    for i, (trid, motifs, alt, gt) in enumerate(rows):
        ref = CAG if motifs != "AT" else "AT" * 8
        pos = 1000 * (i + 1)
        lines.append("\t".join([
            "chr1", str(pos), ".", ref, alt, ".", "PASS",
            f"TRID={trid};END={pos + len(ref) - 1};MOTIFS={motifs}", "GT:AL", f"{gt}:30,30",
        ]))
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    return str(path)


def test_lps_bp():
    assert lps_bp(CAG, "CAG") == 30
    assert lps_bp(INTERRUPTED, "CAG") == 18, "an interruption splits the pure run"
    assert lps_bp("G" + CAG, "AGC") == 30, "any rotation of the motif counts"
    assert lps_bp("TTTT", "CAG") == 0


def test_stats_and_merge(tmp_path: Path):
    # Locus A: LPS per haplotype 30, 18 | 30, 33 | 18, 18.  Locus B: always 16.  Locus C: two motifs, skipped.
    s1 = _vcf(tmp_path / "s1.vcf.gz", "S1", [("A", "CAG", INTERRUPTED, "0/1"), ("B", "AT", ".", "0/0"), ("C", "CAG,CAA", ".", "0/0")])
    s2 = _vcf(tmp_path / "s2.vcf.gz", "S2", [("A", "CAG", CAG + "CAG", "0/1"), ("B", "AT", ".", "0/0"), ("C", "CAG,CAA", ".", "0/0")])
    s3 = _vcf(tmp_path / "s3.vcf.gz", "S3", [("A", "CAG", INTERRUPTED, "1/1"), ("B", "AT", ".", "0/0"), ("C", "CAG,CAA", ".", "./.")])
    cat = tmp_path / "catalog.bed"
    main(["catalog", "--vcf", s1, "--out", str(cat)])
    assert cat.read_text().splitlines()[0] == "chr1\t999\t1029\tA\tCAG"
    sh1, sh2 = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    main(["stats", "--catalog-bed", str(cat), "--vcf", s1, "--vcf", s2, "--out", str(sh1)])
    main(["stats", "--catalog-bed", str(cat), "--vcf", s3, "--threads", "2", "--out", str(sh2)])
    out = tmp_path / "plvi.tsv"
    main(["merge", "--stats", str(sh1), str(sh2), "--out", str(out)])
    rows = {r["trid"]: r for r in csv.DictReader(out.open(), delimiter="\t")}
    assert set(rows) == {"A", "B"}
    vals = [30, 18, 30, 33, 18, 18]
    mean = sum(vals) / 6
    sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / 6)
    assert rows["A"]["n_hap"] == "6" and math.isclose(float(rows["A"]["plvi"]), sd, rel_tol=1e-3)
    assert rows["B"]["plvi"] == "0" and rows["B"]["motif_group"] == "2"
    assert rows["A"]["plvi_pct_in_group"] == "100.0000"
