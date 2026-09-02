#!/usr/bin/env python3
"""Tests for scripts/annotate_repeat_units.py (RU_TEST INFO annotator)."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
FIX = Path(__file__).resolve().parent / "fixtures"

sys.path.insert(0, str(SCRIPTS))

from annotate_repeat_units import (  # noqa: E402
    annotate_record,
    annotate_vcf,
    load_tandem_arrays,
    parse_tandem_bed_line,
)
from sv_site_utils import parse_info as parse_info_shared  # noqa: E402


def _info_map(vcf_text: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for line in vcf_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        vid = parts[2]
        out[vid] = parse_info_shared(parts[7])
    return out


def test_parse_ucsc_and_3col_and_bin():
    three = parse_tandem_bed_line("chr22\t490\t530\n")
    assert three is not None
    chrom, arr = three
    assert chrom == "chr22"
    assert arr.start == 490 and arr.end == 530
    assert arr.period is None and arr.copy_num is None

    ucsc = parse_tandem_bed_line(
        "chr22\t490\t530\ttrf_cag\t3\t10\t3\t90\t2\t300\t20\t30\t40\t10\t1.5\tCAG\n"
    )
    assert ucsc is not None
    _, arr = ucsc
    assert arr.period == 3
    assert arr.copy_num == 10
    assert arr.motif == "CAG"

    dumped = parse_tandem_bed_line("585\tchr22\t7980\t8060\ttrf_a\t1\t40\t1\t95\t0\t100\t100\t0\t0\t0\t0\tA\n")
    assert dumped is not None
    chrom, arr = dumped
    assert chrom == "chr22"
    assert arr.period == 1
    assert arr.copy_num == 40
    assert arr.motif == "A"

    trf = parse_tandem_bed_line("chr22\t100\t200\t4\t12.5\tATGC\n")
    assert trf is not None
    _, arr = trf
    assert arr.period == 4
    assert arr.copy_num == 12.5
    assert arr.motif == "ATGC"


def test_annotate_vcf_marks_only_ru_eligible(tmp_path: Path):
    arrays = load_tandem_arrays(str(FIX / "simpleRepeat.ru.bed"))
    out = tmp_path / "annotated.vcf"
    with (FIX / "ru_sites.vcf").open() as inf, out.open("w") as fh:
        n_ru = annotate_vcf(inf, fh, arrays)
    text = out.read_text()
    assert "ID=RU_TEST" in text
    assert "ID=PERIOD" in text
    assert "ID=CN_REF" in text
    assert "Number=A" in text
    info = _info_map(text)

    assert n_ru == 5
    assert "RU_TEST" in info["vntr_ins"]
    assert info["vntr_ins"]["PERIOD"] == "3"
    assert info["vntr_ins"]["MOTIF"] == "CAG"
    assert info["vntr_ins"]["CN_REF"] == "10"
    assert info["vntr_ins"]["RU"] == "3"

    assert "RU_TEST" in info["vntr_del"]
    assert info["vntr_del"]["PERIOD"] == "2"
    assert info["vntr_del"]["CN_REF"] == "20"
    assert info["vntr_del"]["RU"] == "5"

    assert "RU_TEST" in info["vntr_dup"]
    assert info["vntr_dup"]["PERIOD"] == "5"
    assert info["vntr_dup"]["RU"] == "4"

    assert "RU_TEST" in info["vntr_multi"]
    assert info["vntr_multi"]["RU"] == "1,3"
    assert info["vntr_multi"]["PERIOD"] == "3"

    assert "RU_TEST" in info["hom_ins"]
    assert info["hom_ins"]["PERIOD"] == "1"
    assert info["hom_ins"]["RU"] == "10"
    assert info["hom_ins"]["CN_REF"] == "40"

    for vid in ("snv1", "inv1", "bnd1", "uniq_del", "nohit_ins"):
        assert "RU_TEST" not in info[vid], vid
        assert "PERIOD" not in info[vid], vid


def test_3col_bed_infers_period_no_cn_ref():
    arrays = load_tandem_arrays(str(FIX / "simpleRepeat.ru.3col.bed"))
    ins = annotate_record(
        500, "A", "ACAGCAGCAG", {"SVTYPE": "INS", "SVLEN": "9"}, arrays["chr22"]
    )
    assert ins is not None
    assert ins.period == 3
    assert ins.ru == (3,)
    assert ins.cn_ref is None
    assert ins.motif == "CAG"

    symbolic = annotate_record(
        1000, "N", "<DEL>", {"SVTYPE": "DEL", "SVLEN": "-10", "END": "1010"}, arrays["chr22"]
    )
    assert symbolic is None


def test_inv_bnd_unique_not_ru():
    arrays = load_tandem_arrays(str(FIX / "simpleRepeat.ru.bed"))
    assert (
        annotate_record(3000, "N", "<INV>", {"SVTYPE": "INV", "SVLEN": "140", "END": "3140"}, arrays["chr22"])
        is None
    )
    assert annotate_record(4000, "N", "N]chr1:100]", {"SVTYPE": "BND"}, arrays["chr22"]) is None
    assert (
        annotate_record(5000, "N", "<DEL>", {"SVTYPE": "DEL", "SVLEN": "-17", "END": "5017"}, arrays["chr22"])
        is None
    )
    assert annotate_record(100, "A", "G", {}, arrays.get("chr22", [])) is None


def test_cli_roundtrip(tmp_path: Path):
    import subprocess

    out = tmp_path / "out.vcf.gz"
    subprocess.check_call(
        [
            sys.executable,
            str(SCRIPTS / "annotate_repeat_units.py"),
            "--vcf",
            str(FIX / "ru_sites.vcf"),
            "--simple-repeat-bed",
            str(FIX / "simpleRepeat.ru.bed"),
            "--out",
            str(out),
        ]
    )
    with gzip.open(out, "rt") as fh:
        text = fh.read()
    info = _info_map(text)
    assert "RU_TEST" in info["vntr_ins"]
    assert info["vntr_multi"]["RU"] == "1,3"
