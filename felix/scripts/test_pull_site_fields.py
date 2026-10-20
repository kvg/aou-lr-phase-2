"""Tests for pull_site_fields: the indexed range-read, the site matching, and the two case/control summaries.

The query tests build real bgzipped, indexed VCFs in tmp_path and run real bcftools, so they exercise the thing that
actually talks to htslib. They skip when no bcftools is on PATH (set BCFTOOLS to point at one).
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pull_site_fields as ps  # noqa: E402

BCFTOOLS = os.environ.get("BCFTOOLS") or shutil.which("bcftools")
needs_bcftools = pytest.mark.skipif(not BCFTOOLS, reason="no bcftools on PATH (set BCFTOOLS)")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=RNC,Number=2,Type=Character,Description="Reason for no call">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
"""


def _vcf(tmp_path: Path, rows: list[str], name="x.vcf.gz", header=HEADER) -> str:
    plain = tmp_path / name.replace(".vcf.gz", ".vcf")
    plain.write_text(header + "".join(r if r.endswith("\n") else r + "\n" for r in rows))
    out = tmp_path / name
    subprocess.run([BCFTOOLS, "view", "-Oz", "-o", str(out), str(plain)], check=True,
                   capture_output=True)
    subprocess.run([BCFTOOLS, "index", "-t", str(out)], check=True, capture_output=True)
    os.environ["BCFTOOLS"] = BCFTOOLS
    return str(out)


def _sites(tmp_path: Path, rows: list[dict], name="sites.tsv") -> Path:
    p = tmp_path / name
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


# ------------------------------------------------------------------------------------------- the site list


def test_read_sites_requires_the_three_core_columns(tmp_path):
    p = _sites(tmp_path, [{"locus_id": "a", "chrom": "9"}])
    with pytest.raises(KeyError, match="pos"):
        ps.read_sites(p)


def test_read_sites_fills_absent_allele_columns(tmp_path):
    df = ps.read_sites(_sites(tmp_path, [{"locus_id": "a", "chrom": "9", "pos": 100}]))
    assert df.ref.isna().all() and df.alt.isna().all() and df.pos.dtype.kind == "i"


def test_write_regions_is_one_line_per_site_in_order(tmp_path):
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "b", "chrom": "9", "pos": 200},
                                            {"locus_id": "a", "chrom": "1", "pos": 100}]))
    out = tmp_path / "r.tsv"
    assert ps.write_regions(sites, out, chrom_prefix="chr") == 2
    assert out.read_text().splitlines() == ["chr1\t100\t100", "chr9\t200\t200"]


# -------------------------------------------------------------------------------------------- the querying


@needs_bcftools
def test_available_format_tags_reads_the_header(tmp_path):
    vcf = _vcf(tmp_path, ["chr9\t100\t.\tA\tT\t.\t.\t.\tGT:DP:GQ\t0/0:30:40\t0/1:28:38\t./.:0:0\t1/1:31:41"])
    assert {"GT", "DP", "GQ", "RNC"} <= ps.available_format_tags(vcf)


@needs_bcftools
def test_samples_come_back_in_file_order(tmp_path):
    vcf = _vcf(tmp_path, ["chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1"])
    assert ps.samples(vcf) == ["s1", "s2", "s3", "s4"]


@needs_bcftools
def test_query_sites_returns_one_row_per_site_and_sample(tmp_path):
    vcf = _vcf(tmp_path, [
        "chr9\t100\t.\tA\tT\t.\t.\t.\tGT:DP:GQ\t0/0:30:40\t0/1:28:38\t./.:0:0\t1/1:31:41",
        "chr9\t500\t.\tC\tG\t.\t.\t.\tGT:DP:GQ\t0/0:20:30\t0/0:22:33\t0/1:25:35\t./.:0:0",
        "chr9\t900\t.\tC\tG\t.\t.\t.\tGT:DP:GQ\t0/0:20:30\t0/0:22:33\t0/1:25:35\t0/0:19:29",
    ])
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "L1", "chrom": "9", "pos": 100},
                                            {"locus_id": "L2", "chrom": "9", "pos": 500}]))
    got = ps.query_sites(vcf, sites, fields=("GT", "DP", "GQ"))
    assert len(got) == 8  # two sites requested, four samples; the 900 record is not read
    assert set(got.locus_id) == {"L1", "L2"}
    l1 = got[got.locus_id == "L1"].set_index("sample")
    assert l1.loc["s3", "GT"] == "./." and l1.loc["s4", "GT"] == "1/1"
    assert l1.loc["s1", "DP"] == "30"


@needs_bcftools
def test_query_sites_matches_the_right_record_when_a_position_is_shared(tmp_path):
    """Two records at one POS: a POS-only match has picked the wrong marker in this project before."""
    vcf = _vcf(tmp_path, [
        "chr9\t100\t.\tA\tT\t.\t.\t.\tGT:DP\t0/0:10\t0/0:10\t0/0:10\t0/0:10",
        "chr9\t100\t.\tA\tAT\t.\t.\t.\tGT:DP\t0/1:99\t0/1:99\t0/1:99\t0/1:99",
    ])
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "L1", "chrom": "9", "pos": 100,
                                             "ref": "A", "alt": "AT"}]))
    got = ps.query_sites(vcf, sites, fields=("GT", "DP"))
    assert len(got) == 4 and set(got.DP) == {"99"} and set(got.alt) == {"AT"}


@needs_bcftools
def test_query_sites_without_alleles_takes_every_record_at_the_position(tmp_path):
    vcf = _vcf(tmp_path, [
        "chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0",
        "chr9\t100\t.\tA\tAT\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1",
    ])
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "L1", "chrom": "9", "pos": 100}]))
    assert len(ps.query_sites(vcf, sites, fields=("GT",))) == 8


@needs_bcftools
def test_query_sites_reports_an_undeclared_tag_as_nan_rather_than_failing(tmp_path):
    header = HEADER.replace('##FORMAT=<ID=RNC,Number=2,Type=Character,Description="Reason for no call">\n', "")
    vcf = _vcf(tmp_path, ["chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1"], header=header)
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "L1", "chrom": "9", "pos": 100}]))
    got = ps.query_sites(vcf, sites, fields=("GT", "RNC"))
    assert len(got) == 4 and got.RNC.isna().all() and set(got.GT) == {"0/0", "0/1", "./.", "1/1"}


@needs_bcftools
def test_query_sites_is_empty_when_the_position_holds_nothing(tmp_path):
    vcf = _vcf(tmp_path, ["chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1"])
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "L1", "chrom": "9", "pos": 7_000}]))
    assert ps.query_sites(vcf, sites, fields=("GT",)).empty


@needs_bcftools
def test_query_sites_reads_several_chromosomes_in_one_pass(tmp_path):
    vcf = _vcf(tmp_path, [
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/0\t0/0",
        "chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t0/0\t0/0",
    ])
    sites = ps.read_sites(_sites(tmp_path, [{"locus_id": "A", "chrom": "1", "pos": 100},
                                            {"locus_id": "B", "chrom": "9", "pos": 100}]))
    got = ps.query_sites(vcf, sites, fields=("GT",))
    assert set(got.locus_id) == {"A", "B"} and len(got) == 8


@needs_bcftools
def test_query_sites_raises_with_the_stderr_when_bcftools_fails(tmp_path):
    with pytest.raises(RuntimeError, match="bcftools view -h failed"):
        ps.available_format_tags(str(tmp_path / "nope.vcf.gz"))


# ---------------------------------------------------------------------------------------------- no-calls


@pytest.mark.parametrize("gt,expected", [("./.", True), (".|.", True), (".", True), ("0/0", False),
                                         ("0|1", False), ("1/1", False), (None, True)])
def test_is_nocall_recognises_every_missing_form(gt, expected):
    assert ps.is_nocall(gt) is expected


def _rows(per_sample: dict[str, tuple[str, int, int]], locus="L1"):
    return pd.DataFrame([dict(locus_id=locus, chrom="9", pos=100, ref="A", alt="T", sample=s,
                              GT=g, DP=d, GQ=q) for s, (g, d, q) in per_sample.items()])


def test_call_rate_by_group_counts_each_group_separately():
    rows = _rows({"c1": ("0/0", 30, 40), "c2": ("./.", 0, 0), "n1": ("0/1", 31, 41), "n2": ("0/0", 29, 39)})
    groups = {"c1": "case", "c2": "case", "n1": "control", "n2": "control"}
    out = ps.call_rate_by_group(rows, groups).iloc[0]
    assert out.n_case == 2 and out.n_case_nocall == 1 and out.call_rate_case == 0.5
    assert out.n_control == 2 and out.n_control_nocall == 0 and out.call_rate_control == 1.0
    assert out.call_rate_diff == -0.5
    assert out.median_DP_control == 30.0


def test_call_rate_by_group_ignores_samples_outside_the_phecode():
    rows = _rows({"c1": ("0/0", 30, 40), "x1": ("./.", 0, 0)})
    out = ps.call_rate_by_group(rows, {"c1": "case"}).iloc[0]
    assert out.n_case == 1 and math.isnan(out.call_rate_control) is False or out.n_control == 0


def test_call_rate_by_group_is_empty_without_any_assigned_sample():
    rows = _rows({"c1": ("0/0", 30, 40)})
    assert ps.call_rate_by_group(rows, {"other": "case"}).empty


def test_imputed_fraction_counts_pre_nocalls_that_phasing_filled():
    pre = _rows({"c1": ("./.", 0, 0), "c2": ("0/0", 30, 40), "n1": ("0/0", 30, 40), "n2": ("0/0", 30, 40)})
    post = _rows({"c1": ("0/1", 0, 0), "c2": ("0/0", 30, 40), "n1": ("0/0", 30, 40), "n2": ("0/0", 30, 40)})
    groups = {"c1": "case", "c2": "case", "n1": "control", "n2": "control"}
    out = ps.imputed_fraction_by_group(pre, post, groups).iloc[0]
    assert out.imputed_frac_case == 0.5 and out.imputed_frac_control == 0.0
    assert out.imputed_frac_diff == 0.5 and out.n_compared == 4


def test_imputed_fraction_does_not_count_a_genotype_still_missing_after_phasing():
    pre = _rows({"c1": ("./.", 0, 0)})
    post = _rows({"c1": ("./.", 0, 0)})
    out = ps.imputed_fraction_by_group(pre, post, {"c1": "case"}).iloc[0]
    assert out.imputed_frac_case == 0.0


# ------------------------------------------------------------------------------------- case/control map


def _pheno(tmp_path):
    return pd.DataFrame({"person_id": ["1", "2", "3", "4"], "BI_174": [1, 0, 0, None], "CV_400": [0, 0, 1, 1]})


def test_case_control_map_reads_one_phecode(tmp_path):
    m = ps.case_control_map(_pheno(tmp_path), "BI_174")
    assert m == {"1": "case", "2": "control", "3": "control"}  # person 4 is neither


def test_case_control_map_restricts_to_the_analysis_samples(tmp_path):
    m = ps.case_control_map(_pheno(tmp_path), "BI_174", analysis_samples={"1", "3"})
    assert m == {"1": "case", "3": "control"}


def test_case_control_map_rejects_an_unknown_phecode(tmp_path):
    with pytest.raises(KeyError, match="ZZ_999"):
        ps.case_control_map(_pheno(tmp_path), "ZZ_999")


def test_phenotype_column_prefers_the_first_name_present(tmp_path):
    df = pd.DataFrame({"person_id": ["1"], "BI_164_6": [1]})
    assert ps.phenotype_column(df, ["BI_164_6", "BI_164.6"]) == "BI_164_6"
    assert ps.phenotype_column(df, ["BI_164.6", "BI_164_6"]) == "BI_164_6"   # falls through to the one that exists
    assert ps.phenotype_column(df, "BI_164_6") == "BI_164_6"


def test_phenotype_column_names_everything_it_tried(tmp_path):
    df = pd.DataFrame({"person_id": ["1"], "BI_164_6": [1]})
    with pytest.raises(KeyError, match=r"XX_1.*XX\.1|XX\.1.*XX_1"):
        ps.phenotype_column(df, ["XX_1", "XX.1"])


def test_case_control_map_finds_the_underscore_column_given_both_forms(tmp_path):
    """pheno_cov uses BI_164_6; the figure code and PhecodeX tables use BI_164.6. Both are offered, one is found."""
    df = pd.DataFrame({"person_id": ["1", "2", "3"], "BI_164_6": [1, 0, None]})
    assert ps.case_control_map(df, ["BI_164_6", "BI_164.6"]) == {"1": "case", "2": "control"}
    assert ps.case_control_map(df, ["BI_164.6", "BI_164_6"]) == {"1": "case", "2": "control"}


def test_imputed_fraction_without_the_phased_rows_counts_every_pre_nocall():
    """Phasing leaves no no-calls behind, so the phased query is unnecessary; passing None asserts that."""
    pre = _rows({"c1": ("./.", 0, 0), "c2": ("0/0", 30, 40), "n1": ("0/0", 30, 40), "n2": ("0/0", 30, 40)})
    groups = {"c1": "case", "c2": "case", "n1": "control", "n2": "control"}
    out = ps.imputed_fraction_by_group(pre, None, groups).iloc[0]
    assert out.imputed_frac_case == 0.5 and out.imputed_frac_control == 0.0


def test_imputed_fraction_with_phased_rows_can_contradict_the_assumption():
    """If a genotype is still a no-call after phasing it is not an imputed genotype, and passing post shows that."""
    pre = _rows({"c1": ("./.", 0, 0)})
    post = _rows({"c1": ("./.", 0, 0)})
    assert ps.imputed_fraction_by_group(pre, post, {"c1": "case"}).iloc[0].imputed_frac_case == 0.0
    assert ps.imputed_fraction_by_group(pre, None, {"c1": "case"}).iloc[0].imputed_frac_case == 1.0


@needs_bcftools
def test_probe_reports_success_on_a_readable_file(tmp_path):
    vcf = _vcf(tmp_path, ["chr9\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1"])
    assert ps.probe(vcf).startswith("OK:")


@needs_bcftools
def test_probe_reports_what_bcftools_said_on_failure(tmp_path):
    out = ps.probe(str(tmp_path / "absent.vcf.gz"))
    assert out.startswith("FAILED") and "bcftools said" in out
    assert "environment carrying credentials" in out
