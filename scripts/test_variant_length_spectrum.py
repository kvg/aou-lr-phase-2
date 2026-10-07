#!/usr/bin/env python3
"""Tiny synthetic GLnexus shard + SV companions → 1 bp × region_class counts."""

from __future__ import annotations

import gzip
import json
import subprocess
import tempfile
from pathlib import Path

from annotate_repeat_context import context_label
from summarize_variant_length_spectrum import (
    BINS_HEADER,
    breakpoint_class,
    chrom_contig,
    chrom_matches,
    small_allele_length,
)

SCRIPTS = Path(__file__).resolve().parent
HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr22>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def test_helpers() -> None:
    assert chrom_contig("aou_lr_phase2_v1.chr1") == "chr1"
    assert chrom_contig("chr1") == "chr1"
    assert chrom_matches("1", "chr1")
    assert not chrom_matches("chr2", "chr1")
    assert small_allele_length("A", "G") == ("snv", 0)
    assert small_allele_length("AT", "GT") == ("snv", 0)
    assert small_allele_length("AT", "GC") == ("mnp", None)
    assert small_allele_length("A", "*") == ("other", None)
    assert small_allele_length("A", "ATT") == ("indel", 2)
    assert small_allele_length("ATT", "A") == ("indel", -2)


def test_breakpoint_class_matches_context_label() -> None:
    iv = {"RM": [(10, 20), (40, 60)], "SD": [(50, 80)], "SR": [(55, 58)]}
    points = [([s for s, _ in iv[k]], [e for _, e in iv[k]]) for k in ("RM", "SD", "SR")]
    for start0 in range(0, 90):
        for width in (1, 2, 5, 30):
            want, _ = context_label("SNV", start0, start0 + width, iv)
            assert breakpoint_class(points, start0, start0 + width) == want, (start0, width)


def test_end_to_end_script() -> None:
    glnexus = HEADER + "".join(
        [
            "chr22\t100\t.\tA\tG\t.\t.\t.\n",            # SNV, US
            "chr22\t1001\t.\tA\tG,AT\t.\tPASS\t.\n",     # SNV + INS1 in RM
            "chr22\t1101\t.\tAT\tGC\t.\t.\t.\n",         # MNP: summary only
            "chr22\t1201\t.\tA\t*,C\t.\t.\t.\n",         # * skipped, SNV kept
            "chr22\t300\t.\tA\tT\t.\tLowQual\t.\n",      # filtered
            "chr22\t400\t.\t" + "A" * 20 + "\tA\t.\t.\t.\n",  # DEL19 kept
            "chr22\t500\t.\t" + "A" * 21 + "\tA\t.\t.\t.\n",  # DEL20 dropped (SV source)
            "chr1\t10\t.\tA\tG\t.\t.\t.\n",              # off-contig
        ]
    )
    main = HEADER + "".join(
        [
            "chr22\t2000\tdel_alu\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-310;END=2310\n",
            "chr22\t3000\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=6000\n",
            "chr22\t4000\tsmall\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-19;END=4019\n",
            "chr22\t4100\tins20\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=20\n",
            "chr22\t5000\tinv\tN\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=500;END=5500\n",
            "chr1\t5000\tother\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=80\n",
        ]
    )
    large = HEADER + "chr22\t7000\tul\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-25000;END=32000\n"
    bnd = HEADER + "".join(
        [
            "chr22\t600\tbnd1\tN\tN]chr22:1000]\t.\tPASS\tSVTYPE=BND;MATEID=bnd2\n",  # mate end in RM
            "chr22\t1000\tbnd2\tN\t]chr22:600]N\t.\tPASS\tSVTYPE=BND;MATEID=bnd1\n",  # own end in RM
            "chr22\t8000\tbnd3\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr1;END=100\n",   # own end in SD
            "chr22\t600\tbnd4\tN\t<BND>\t.\tPASS\tSVTYPE=BND\n",                      # no mate: US
        ]
    )
    beds = {
        "rmsk": "chr22\t990\t1010\n",
        "sr": "chr22\t2305\t2320\n",
        "sd": "chr22\t7000\t40000\n",
    }
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        for name, text in (("gl", glnexus), ("main", main), ("large", large), ("bnd", bnd)):
            (t / f"{name}.vcf").write_text(text, encoding="utf-8")
        for name, text in beds.items():
            with gzip.open(t / f"{name}.bed.gz", "wt") as fh:
                fh.write(text)
        bed_args = [
            "--rmsk-bed", str(t / "rmsk.bed.gz"),
            "--simple-repeat-bed", str(t / "sr.bed.gz"),
            "--segdup-bed", str(t / "sd.bed.gz"),
        ]
        run = ["python3", str(SCRIPTS / "summarize_variant_length_spectrum.py")]
        subprocess.check_call(
            run + ["--glnexus-vcf", str(t / "gl.vcf"), "--chrom", "chr22", *bed_args,
                   "--out-bins", str(t / "gl.tsv"), "--out-summary", str(t / "gl.json")]
        )
        subprocess.check_call(
            run + ["--main-vcf", str(t / "main.vcf"), "--large-vcf", str(t / "large.vcf"),
                   "--bnd-vcf", str(t / "bnd.vcf"), *bed_args,
                   "--out-bins", str(t / "comp.tsv"), "--out-summary", str(t / "comp.json")]
        )
        subprocess.check_call(
            ["python3", str(SCRIPTS / "merge_variant_length_spectrum.py"),
             "--bins", str(t / "gl.tsv"), str(t / "comp.tsv"),
             "--summaries", str(t / "gl.json"), str(t / "comp.json"),
             "--out-bins", str(t / "all.tsv"), "--out-summary", str(t / "all.json")]
        )

        gl = (t / "gl.tsv").read_text(encoding="utf-8")
        assert gl.startswith(BINS_HEADER)
        assert "chr22\tsmall\t0\tUS\t2\n" in gl
        assert "chr22\tsmall\t0\tRM\t1\n" in gl
        assert "chr22\tsmall\t1\tRM\t1\n" in gl
        assert "chr22\tsmall\t-19\tUS\t1\n" in gl
        assert "\t-20\t" not in gl
        g = json.loads((t / "gl.json").read_text(encoding="utf-8"))["counts"]
        assert g["snv_records"] == 3
        assert g["glnexus_records_filtered"] == 1
        assert g["glnexus_records_off_contig"] == 1
        assert g["glnexus_alleles_mnp"] == 1
        assert g["glnexus_alleles_other"] == 1
        assert g["glnexus_indel_ge_small_max_dropped"] == 1

        comp = (t / "comp.tsv").read_text(encoding="utf-8")
        assert "companions\tsv\t-310\tSR\t1\n" in comp
        assert "companions\tsv\t6000\tUS\t1\n" in comp
        assert "companions\tsv\t20\tUS\t1\n" in comp
        assert "companions\tsv\t80\tUS\t1\n" in comp
        assert "\t-19\t" not in comp
        assert "\t500\t" not in comp
        assert "companions\tultralong\t-25000\tSD\t1\n" in comp
        c = json.loads((t / "comp.json").read_text(encoding="utf-8"))["counts"]
        assert c["main_ge20_DEL"] == 1 and c["main_ge50_DEL"] == 1
        assert c["main_ge20_INS"] == 3 and c["main_ge50_INS"] == 2
        assert c["main_ge20_INV"] == 1
        assert c["main_lt_sv_min_dropped"] == 1
        assert c["bnd_records"] == 4 and c["large_records"] == 1
        assert c["bnd_mate_unparsed"] == 1
        assert c["bnd_records_with_mateid"] == 2 and c["bnd_records_mate_in_file"] == 2
        assert "companions\tbnd\t0\tRM\t2\n" in comp
        assert "companions\tbnd\t0\tSD\t1\n" in comp
        assert "companions\tbnd\t0\tUS\t1\n" in comp

        merged = json.loads((t / "all.json").read_text(encoding="utf-8"))
        assert merged["table2_check"]["SNVs (bcftools 'number of SNPs')"] == 3
        assert merged["table2_check"]["BND"] == 4
        assert merged["partition_totals"] == {"bnd": 4, "small": 5, "sv": 4, "ultralong": 1}


if __name__ == "__main__":
    test_helpers()
    test_breakpoint_class_matches_context_label()
    test_end_to_end_script()
    print("ok")
