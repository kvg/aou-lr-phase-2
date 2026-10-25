"""Tests for trgt_polymorphism_spectrum.py."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import trgt_polymorphism_spectrum as tps  # noqa: E402


def fold(al_values, n_loci=None):
    """Fold one sample's AL column into a fresh accumulator."""
    series = pd.Series(al_values, dtype=str)
    acc = tps.new_accumulator(n_loci or len(series))
    tps.accumulate_frame(acc, 0, series)
    return acc


def test_heterozygous_call_counted_once():
    acc = fold(["100,130"])
    assert acc["n_called"][0] == 1
    assert acc["n_alleles"][0] == 2
    assert acc["n_het"][0] == 1
    assert acc["al_min"][0] == 100 and acc["al_max"][0] == 130
    assert acc["al_sum"][0] == 230
    assert acc["al_sumsq"][0] == 100 ** 2 + 130 ** 2


def test_equal_lengths_are_not_heterozygous():
    """The same length twice is homozygous even though the GT may be 0/1.

    This is the whole point of the metric: a dosage test needs the two
    haplotypes to differ in LENGTH, not merely in allele identity.
    """
    acc = fold(["115,115"])
    assert acc["n_het"][0] == 0
    assert acc["n_alleles"][0] == 2
    assert acc["al_min"][0] == acc["al_max"][0] == 115


def test_hemizygous_call_has_one_allele_and_no_het():
    acc = fold(["97"])
    assert acc["n_called"][0] == 1
    assert acc["n_alleles"][0] == 1
    assert acc["n_het"][0] == 0
    assert acc["al_min"][0] == acc["al_max"][0] == 97


@pytest.mark.parametrize("missing", [".", "", ".,."])
def test_no_call_contributes_nothing(missing):
    acc = fold([missing])
    assert acc["n_called"][0] == 0
    assert acc["n_alleles"][0] == 0
    assert acc["n_het"][0] == 0
    # a locus never called must not look like a zero-length locus
    assert acc["al_min"][0] == tps.MISSING_MAX
    assert acc["al_max"][0] == -1


def test_half_missing_call_keeps_the_observed_allele():
    acc = fold(["120,."])
    assert acc["n_called"][0] == 1
    assert acc["n_alleles"][0] == 1
    assert acc["n_het"][0] == 0, "one observed allele cannot establish heterozygosity"
    assert acc["al_min"][0] == acc["al_max"][0] == 120


def test_merge_is_equivalent_to_sequential_folding():
    seq = tps.new_accumulator(2)
    for values in (["100,130", "50,50"], ["140,100", "50,55"]):
        tps.accumulate_frame(seq, 0, pd.Series(values, dtype=str))
    merged = tps.merge(fold(["100,130", "50,50"]), fold(["140,100", "50,55"]))
    for field in tps.FIELDS:
        assert np.array_equal(seq[field], merged[field]), field
    assert seq["al_min"][0] == 100 and seq["al_max"][0] == 140
    assert seq["n_het"][1] == 1 and seq["n_called"][1] == 2


def test_offset_places_chunks_at_their_own_loci():
    acc = tps.new_accumulator(4)
    tps.accumulate_frame(acc, 0, pd.Series(["10,10"], dtype=str))
    tps.accumulate_frame(acc, 2, pd.Series(["20,30"], dtype=str))
    assert acc["n_called"].tolist() == [1, 0, 1, 0]
    assert acc["n_het"].tolist() == [0, 0, 1, 0]


def test_summary_separates_monomorphic_from_uncalled():
    acc = tps.new_accumulator(3)
    # locus 0 polymorphic, locus 1 monomorphic, locus 2 never called
    for values in (["100,130", "50,50", "."], ["110,100", "50,50", "."]):
        tps.accumulate_frame(acc, 0, pd.Series(values, dtype=str))
    summary = tps.summarise(acc, [0.5], cohort_n=12261).set_index("metric")["value"]
    assert summary["loci_total"] == 3
    assert summary["loci_with_any_call"] == 2
    assert summary["loci_monomorphic"] == 1
    assert summary["loci_polymorphic"] == 1
    assert summary["loci_het_rate_ge_0.5"] == 1


def test_uncalled_locus_is_never_reported_as_polymorphic():
    acc = tps.new_accumulator(1)
    tps.accumulate_frame(acc, 0, pd.Series(["."], dtype=str))
    summary = tps.summarise(acc, [], cohort_n=1).set_index("metric")["value"]
    assert summary["loci_polymorphic"] == 0
    assert summary["loci_monomorphic"] == 0


def _fake_bcftools(tmp_path: Path, body: str, rc: int = 0) -> str:
    script = tmp_path / "fake_bcftools"
    script.write_text("#!/bin/sh\ncat <<'EOF'\n" + body + "EOF\nexit %d\n" % rc)
    script.chmod(0o755)
    return str(script)


def test_read_one_hashes_trid_and_folds(tmp_path):
    body = "chr1_100_A\t100,130\nchr1_200_B\t50,50\n"
    acc = tps.new_accumulator(2)
    digest = tps.read_one("x.vcf.gz", acc, 2, bcftools=_fake_bcftools(tmp_path, body))
    assert len(digest) == 32
    assert acc["n_het"].tolist() == [1, 0]


def test_read_one_rejects_a_short_file(tmp_path):
    acc = tps.new_accumulator(5)
    with pytest.raises(SystemExit, match="different catalog or a truncated object"):
        tps.read_one("x.vcf.gz", acc, 5,
                     bcftools=_fake_bcftools(tmp_path, "chr1_100_A\t100,130\n"))


def test_read_one_rejects_a_long_file(tmp_path):
    body = "".join(f"chr1_{i}_A\t10,10\n" for i in range(3))
    acc = tps.new_accumulator(2)
    with pytest.raises(SystemExit, match="different catalog"):
        tps.read_one("x.vcf.gz", acc, 2, bcftools=_fake_bcftools(tmp_path, body))


def test_read_one_propagates_bcftools_failure(tmp_path):
    acc = tps.new_accumulator(1)
    with pytest.raises(SystemExit, match="exited 3"):
        tps.read_one("x.vcf.gz", acc, 1,
                     bcftools=_fake_bcftools(tmp_path, "chr1_100_A\t10,10\n", rc=3))


def test_differing_trid_order_is_rejected(tmp_path, monkeypatch, capsys):
    """Two files whose loci differ must abort, not accumulate by ordinal."""
    a = tmp_path / "a.txt"
    a.write_text("one.vcf.gz\ntwo.vcf.gz\n")
    bodies = {"one.vcf.gz": "chr1_100_A\t100,130\n", "two.vcf.gz": "chr1_999_Z\t100,130\n"}
    real = tps.read_one

    def fake(vcf, acc, n_loci, bcftools="bcftools", chunk_rows=500_000,
             gcs_cat="gcloud storage cat"):
        return real(vcf, acc, n_loci, bcftools=_fake_bcftools(tmp_path, bodies[vcf]))

    monkeypatch.setattr(tps, "read_one", fake)
    monkeypatch.setattr(sys, "argv",
                        ["x", "--vcf-list", str(a), "--n-loci", "1", "--max-samples", "2",
                         "--threads", "1"])
    with pytest.raises(SystemExit, match="disagree about the catalog"):
        tps.main()


def test_sampling_is_reproducible_for_a_seed(tmp_path):
    import random
    paths = [f"s{i}.vcf.gz" for i in range(50)]
    first = random.Random(7).sample(paths, 5)
    assert first == random.Random(7).sample(paths, 5)
    assert first != random.Random(8).sample(paths, 5)


def test_end_to_end_summary(tmp_path, monkeypatch):
    listing = tmp_path / "list.txt"
    listing.write_text("one.vcf.gz\ntwo.vcf.gz\n")
    rows = {
        "one.vcf.gz": "L1\t100,130\nL2\t50,50\nL3\t.\n",
        "two.vcf.gz": "L1\t100,100\nL2\t50,50\nL3\t.\n",
    }
    real = tps.read_one

    def fake(vcf, acc, n_loci, bcftools="bcftools", chunk_rows=500_000,
             gcs_cat="gcloud storage cat"):
        # same TRID column in both, so the digests agree
        return real(vcf, acc, n_loci, bcftools=_fake_bcftools(tmp_path, rows[vcf]))

    monkeypatch.setattr(tps, "read_one", fake)
    out = tmp_path / "summary.tsv"
    monkeypatch.setattr(sys, "argv", ["x", "--vcf-list", str(listing), "--n-loci", "3",
                                      "--thresholds", "0.5", "--out-summary", str(out),
                                      "--threads", "1"])
    assert tps.main() == 0
    got = pd.read_csv(out, sep="\t").set_index("metric")["value"]
    assert got["loci_with_any_call"] == 2
    assert got["loci_monomorphic"] == 1     # L2
    assert got["loci_polymorphic"] == 1     # L1, het in one of two samples
    assert got["loci_het_rate_ge_0.5"] == 1


def test_cli_runs_as_a_script():
    proc = subprocess.run([sys.executable,
                           str(Path(__file__).resolve().parent / "trgt_polymorphism_spectrum.py"),
                           "--help"], capture_output=True, text=True)
    assert proc.returncode == 0
    assert "--vcf-list" in proc.stdout


def _fake_cat(tmp_path: Path, payload: Path, rc: int = 0) -> str:
    """Stand-in for `gcloud storage cat`: ignores the gs:// arg, emits a file."""
    script = tmp_path / "fake_cat"
    script.write_text(f"#!/bin/sh\ncat '{payload}'\nexit {rc}\n")
    script.chmod(0o755)
    return str(script)


def _real_trgt_vcf(tmp_path: Path) -> Path:
    """A TRGT-shaped bgzipped VCF, so real bcftools parses a real file."""
    body = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
chr1\t100\t.\tACG\tA\t0\t.\tTRID=chr1_100_ACG\tGT:AL\t0/1:100,130
chr1\t200\t.\tAT\t.\t0\t.\tTRID=chr1_200_AT\tGT:AL\t0/0:50,50
chr1\t300\t.\tGGC\t.\t0\t.\tTRID=chr1_300_GGC\tGT:AL\t./.:.
"""
    plain = tmp_path / "t.vcf"
    plain.write_text(body)
    subprocess.run(["bgzip", "-kf", str(plain)], check=True)
    return tmp_path / "t.vcf.gz"


@pytest.mark.skipif(not shutil.which("bcftools") or not shutil.which("bgzip"),
                    reason="needs real bcftools/bgzip")
def test_real_bcftools_local_file(tmp_path):
    """Guards the query format against bcftools itself, not against a fake."""
    acc = tps.new_accumulator(3)
    tps.read_one(str(_real_trgt_vcf(tmp_path)), acc, 3)
    assert acc["n_called"].tolist() == [1, 1, 0]
    assert acc["n_het"].tolist() == [1, 0, 0]
    assert acc["al_min"][0] == 100 and acc["al_max"][0] == 130


@pytest.mark.skipif(not shutil.which("bcftools") or not shutil.which("bgzip"),
                    reason="needs real bcftools/bgzip")
def test_gs_path_is_streamed_not_opened(tmp_path):
    """A gs:// path must go through the cat command; htslib cannot open it."""
    vcf = _real_trgt_vcf(tmp_path)
    acc = tps.new_accumulator(3)
    tps.read_one("gs://bucket/t.vcf.gz", acc, 3,
                 gcs_cat=_fake_cat(tmp_path, vcf))
    assert acc["n_called"].tolist() == [1, 1, 0]
    assert acc["n_het"].tolist() == [1, 0, 0]


@pytest.mark.skipif(not shutil.which("bcftools"), reason="needs real bcftools")
def test_failed_stream_is_not_silently_short(tmp_path):
    """A broken download must raise, not report fewer loci."""
    vcf = _real_trgt_vcf(tmp_path)
    acc = tps.new_accumulator(3)
    with pytest.raises(SystemExit, match="streaming the object exited|different catalog"):
        tps.read_one("gs://bucket/t.vcf.gz", acc, 3,
                     gcs_cat=_fake_cat(tmp_path, vcf, rc=7))


def test_open_query_uses_cat_only_for_remote_paths(tmp_path, monkeypatch):
    seen = []

    class FakeStream:
        def close(self):
            pass

    class FakePopen:
        def __init__(self, cmd, **kw):
            seen.append(cmd)
            self.stdout, self.stderr = FakeStream(), FakeStream()

    monkeypatch.setattr(tps.subprocess, "Popen", FakePopen)
    tps.open_query("/local/x.vcf.gz", "bcftools", "gcloud storage cat")
    assert len(seen) == 1 and seen[0][0] == "bcftools" and seen[0][-1] == "/local/x.vcf.gz"
    seen.clear()
    tps.open_query("gs://b/x.vcf.gz", "bcftools", "gcloud storage cat")
    assert len(seen) == 2
    assert seen[0] == ["gcloud", "storage", "cat", "gs://b/x.vcf.gz"]
    assert seen[1][-1] == "-", "bcftools must read the stream from stdin"


@pytest.mark.skipif(not shutil.which("bcftools") or not shutil.which("bgzip"),
                    reason="needs real bcftools/bgzip")
def test_read_all_serial_matches_sequential_folding(tmp_path):
    vcf = str(_real_trgt_vcf(tmp_path))
    total = tps.read_all([vcf, vcf, vcf], 3, "bcftools", "gcloud storage cat", threads=1)
    assert total["n_called"].tolist() == [3, 3, 0]
    assert total["n_het"].tolist() == [3, 0, 0]
    assert total["al_sum"][0] == 3 * (100 + 130)


def test_read_all_rejects_mixed_catalogs(tmp_path, monkeypatch):
    """Two layouts must abort even though each file is internally consistent."""
    layouts = iter(["aaaa", "bbbb"])

    def fake(vcf, acc, n_loci, bcftools="bcftools", gcs_cat="x"):
        return next(layouts)

    monkeypatch.setattr(tps, "read_one", fake)
    with pytest.raises(SystemExit, match="do not share a TRID column"):
        tps.read_all(["a.vcf.gz", "b.vcf.gz"], 1, "bcftools", "x", threads=1)
