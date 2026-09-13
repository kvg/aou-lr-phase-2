#!/usr/bin/env python3
"""Unit tests for flare_split_samples and flare_summarize_models (stdlib only)."""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from flare_split_samples import (  # noqa: E402
    load_covariate_rows,
    main as split_main,
    split_samples,
)
from flare_resolve_region import match_region_chrom, parse_contigs_from_header  # noqa: E402
from flare_summarize_models import (  # noqa: E402
    merge_global_anc,
    parse_model,
    pop_from_model_path,
    write_models_tsv,
)


def _write_samples(path: Path, samples: list[str]) -> None:
    path.write_text("".join(f"{s}\n" for s in samples))


def test_split_keeps_pops_and_drops_small():
    vcf = ["1", "2", "3", "4", "5", "6", "7"]
    cov_rows = {
        "1": ("AFR", False),
        "2": ("AFR", False),
        "3": ("AFR", False),
        "4": ("EUR", False),
        "5": ("EUR", False),
        "6": ("MID", False),
        "7": (None, False),
    }
    result = split_samples(
        vcf,
        cov_rows,
        drop_pops=["MID"],
        min_samples=3,
    )
    assert set(result["kept"]) == {"AFR"}
    assert result["kept"]["AFR"] == ["1", "2", "3"]
    assert result["skipped_small"] == {"EUR": 2}
    assert result["unmatched_vcf"] == ["7"]
    assert result["dropped"] == [("6", "MID", "drop_pop")]


def test_include_pops_and_unlabeled():
    result = split_samples(
        ["a", "b", "c", "d"],
        {"a": ("afr", False), "b": ("AMR", False), "c": (None, False), "d": ("EUR", False)},
        include_pops=["AFR", "OTH"],
        unlabeled_pop="OTH",
        min_samples=1,
    )
    assert result["kept"]["AFR"] == ["a"]
    assert result["kept"]["OTH"] == ["c"]
    assert "EUR" not in result["kept"]


def test_excludes_controls_by_flag_and_hg_id():
    result = split_samples(
        ["1", "2", "HG00096", "3"],
        {
            "1": ("AFR", False),
            "2": ("AFR", True),
            "HG00096": ("EUR", False),
            "3": ("AFR", False),
        },
        min_samples=1,
        exclude_controls=True,
    )
    assert result["kept"]["AFR"] == ["1", "3"]
    assert result["excluded_controls"] == ["2", "HG00096"]
    kept = split_samples(
        ["1", "2", "HG00096"],
        {"1": ("AFR", False), "2": ("AFR", True), "HG00096": ("EUR", False)},
        min_samples=1,
        exclude_controls=False,
    )
    assert kept["kept"]["AFR"] == ["1", "2"]
    assert kept["kept"]["EUR"] == ["HG00096"]


def test_cli_gzip_covariates_and_fail_unmatched():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        cov = tmp / "cov.csv.gz"
        with gzip.open(cov, "wt") as fh:
            fh.write(
                "research_id,population,is_reference_control\n"
                "1,AFR,False\n2,AFR,False\nHG00096,EUR,True\n"
            )
        samples = tmp / "samples.txt"
        _write_samples(samples, ["1", "2", "3", "HG00096"])
        out = tmp / "out"
        rc = split_main(
            [
                "--covariates",
                str(cov),
                "--vcf-samples",
                str(samples),
                "--out-dir",
                str(out),
                "--min-samples",
                "1",
                "--max-unmatched-vcf-frac",
                "0.5",
            ]
        )
        assert rc == 0
        assert (out / "keeps" / "AFR.samples.txt").read_text() == "1\n2\n"
        assert (out / "excluded_controls.txt").read_text() == "HG00096\n"
        try:
            split_main(
                [
                    "--covariates",
                    str(cov),
                    "--vcf-samples",
                    str(samples),
                    "--out-dir",
                    str(tmp / "out2"),
                    "--min-samples",
                    "1",
                    "--max-unmatched-vcf-frac",
                    "0.01",
                ]
            )
        except SystemExit as exc:
            assert "unmatched" in str(exc)
        else:
            raise AssertionError("expected unmatched-fraction failure")


def test_load_covariate_rows_duplicate():
    with tempfile.TemporaryDirectory() as td:
        cov = Path(td) / "cov.csv"
        cov.write_text("research_id,population,age\n1,AFR,40\n1,EUR,40\n")
        try:
            load_covariate_rows(cov, "research_id", "population", "is_reference_control")
        except SystemExit as exc:
            assert "duplicate" in str(exc)
        else:
            raise AssertionError("expected duplicate ID failure")


def test_parse_model_and_merge_global():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        model = tmp / "aou.chr22.AFR.model"
        # Full 2A+5 layout (A=5, P=5 → 15 data lines).
        model.write_text(
            "# list of ancestries\n"
            "eas\tamr\teur\tafr\tsas\n"
            "# list of reference panels\n"
            "eas\tamr\teur\tafr\tsas\n"
            "# T\n"
            "11.4\n"
            "# admixture proportions\n"
            "0.01\t0.05\t0.20\t0.70\t0.04\n"
            "1\t0\t0\t0\t0\n"
            "0\t1\t0\t0\t0\n"
            "0\t0\t1\t0\t0\n"
            "0\t0\t0\t1\t0\n"
            "0\t0\t0\t0\t1\n"
            "0.001\t0.01\t0.01\t0.01\t0.01\n"
            "0.01\t0.001\t0.01\t0.01\t0.01\n"
            "0.01\t0.01\t0.001\t0.01\t0.01\n"
            "0.01\t0.01\t0.01\t0.001\t0.01\n"
            "0.01\t0.01\t0.01\t0.01\t0.001\n"
            "0.5\t0.5\t0.5\t0.5\t0.5\n"
        )
        parsed = parse_model(model)
        assert abs(parsed.t_gen - 11.4) < 1e-9
        assert abs(parsed.props[3] - 0.70) < 1e-9
        assert pop_from_model_path(model) == "AFR"

        tsv = tmp / "models.tsv"
        write_models_tsv([("AFR", parsed)], tsv)
        header, row = tsv.read_text().splitlines()
        assert header.startswith("population\tt_gen\tprop_eas")
        assert "mu_eas" in header
        assert row.startswith("AFR\t11.4\t")

        g1 = tmp / "AFR.global.anc.gz"
        g2 = tmp / "EUR.global.anc.gz"
        with gzip.open(g1, "wt") as fh:
            fh.write("SAMPLE eas amr eur afr sas\nS1 0 0 0.1 0.9 0\n")
        with gzip.open(g2, "wt") as fh:
            fh.write("SAMPLE eas amr eur afr sas\nS2 0 0 1 0 0\n")
        merged = tmp / "merged.global.anc.gz"
        n = merge_global_anc([g1, g2], merged)
        assert n == 2
        with gzip.open(merged, "rt") as fh:
            body = fh.read().splitlines()
        assert body[0].startswith("SAMPLE")
        assert body[1].startswith("S1")
        assert body[2].startswith("S2")


def test_match_region_chrom_aliases():
    assert match_region_chrom("chr22:1-10", ["chr22"])[0] == "chr22:1-10"
    assert match_region_chrom("22:1-10", ["chr22"])[0] == "chr22:1-10"
    assert match_region_chrom("chr22:1-10", ["22"])[0] == "22:1-10"
    header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n"
    assert parse_contigs_from_header(header) == ["chr1"]
    try:
        match_region_chrom("chr22:1-10", ["chr1"])
    except SystemExit as exc:
        assert "chr22" in str(exc)
        assert "chr1" in str(exc)
    else:
        raise AssertionError("expected chrom mismatch")


if __name__ == "__main__":
    test_split_keeps_pops_and_drops_small()
    test_include_pops_and_unlabeled()
    test_excludes_controls_by_flag_and_hg_id()
    test_cli_gzip_covariates_and_fail_unmatched()
    test_load_covariate_rows_duplicate()
    test_parse_model_and_merge_global()
    test_match_region_chrom_aliases()
    print("ok")
