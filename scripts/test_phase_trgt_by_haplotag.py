#!/usr/bin/env python3
"""Tests for phase_trgt_by_haplotag.py.

The fixtures are a hand-written TRGT-shaped VCF and a SAM whose reads carry an
allele tag, so the real samtools/bcftools paths are exercised rather than
mocked. Each locus in the fixture is one of the resolution outcomes the script
has to distinguish: oriented from reads, homozygous, conflicting, too thin.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_trgt_by_haplotag as ptb  # noqa: E402

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=TRID,Number=1,Type=String,Description="Locus id">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Allele length">
##FORMAT=<ID=SD,Number=.,Type=Integer,Description="Spanning depth">
##FORMAT=<ID=MC,Number=.,Type=String,Description="Motif counts">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""

# POS, END, TRID, motif, GT, AL, SD, MC
ROWS = [
    # oriented from reads, and the reads say the slots are SWAPPED
    (1000, 1100, "L_swap", "AT", "0/1", "101,151", "5,4", "50,75"),
    # oriented from reads, slots already in haplotype order
    (3000, 3100, "L_keep", "AT", "0/1", "101,151", "5,4", "50,75"),
    # homozygous: nothing to orient
    (5000, 5100, "L_hom", "AT", "0/0", "101,101", "5,5", "50,50"),
    # reads disagree: both slots vote the same haplotype
    (7000, 7100, "L_conflict", "AT", "0/1", "101,151", "5,4", "50,75"),
    # only one haplotagged read per allele
    (9000, 9100, "L_thin", "AT", "0/1", "101,151", "1,1", "50,75"),
    # two ALT alleles: slots carry 1 and 2, reads still carry slot 0 and 1
    (11000, 11100, "L_multi", "AT", "1/2", "121,181", "5,4", "60,90"),
]

# read name -> (locus POS, allele index, haplotype in the haplotag list)
READS = [
    ("r_swap_a1", 1000, 0, 2), ("r_swap_a2", 1000, 0, 2),
    ("r_swap_b1", 1000, 1, 1), ("r_swap_b2", 1000, 1, 1),
    ("r_keep_a1", 3000, 0, 1), ("r_keep_a2", 3000, 0, 1),
    ("r_keep_b1", 3000, 1, 2), ("r_keep_b2", 3000, 1, 2),
    ("r_hom_a1", 5000, 0, 1), ("r_hom_a2", 5000, 0, 2),
    ("r_conf_a1", 7000, 0, 1), ("r_conf_a2", 7000, 0, 1),
    ("r_conf_b1", 7000, 1, 1), ("r_conf_b2", 7000, 1, 1),
    ("r_thin_a1", 9000, 0, 1),
    ("r_thin_b1", 9000, 1, 2),
    ("r_multi_a1", 11000, 0, 2), ("r_multi_a2", 11000, 0, 2),
    ("r_multi_b1", 11000, 1, 1), ("r_multi_b2", 11000, 1, 1),
]


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    vcf = tmp_path / "trgt.vcf"
    lines = [VCF_HEADER.rstrip("\n")]
    for pos, end, trid, motif, gt, al, sd, mc in ROWS:
        info = f"TRID={trid};END={end};MOTIFS={motif}"
        lines.append("\t".join((
            "chr1", str(pos), ".", "A" * 10, "<TR>", ".", ".", info, "GT:AL:SD:MC",
            f"{gt}:{al}:{sd}:{mc}")))
    vcf.write_text("\n".join(lines) + "\n")

    sam = tmp_path / "spanning.sam"
    trid_of = {pos: trid for pos, _e, trid, _m, _g, _a, _s, _c in ROWS}
    sam_lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:248956422"]
    for name, pos, allele, _hap in READS:
        seq = "A" * 200
        sam_lines.append("\t".join((
            name, "0", "chr1", str(pos - 10), "60", "200M", "*", "0", "0", seq,
            "*" * 200, f"TR:Z:{trid_of[pos]}", f"AL:i:{allele}", "FL:i:5")))
    sam.write_text("\n".join(sam_lines) + "\n")

    ht = tmp_path / "haplotag.tsv"
    ht_lines = ["#readname\thaplotype\tphaseset\tchromosome"]
    for name, _pos, _allele, hap in READS:
        ht_lines.append(f"{name}\tH{hap}\t1000\tchr1")
    ht_lines.append("r_untagged\tnone\t.\tchr1")
    ht.write_text("\n".join(ht_lines) + "\n")
    return tmp_path


def run_script(fixture_dir: Path, *extra: str) -> dict:
    out_summary = fixture_dir / "summary.json"
    cmd = [sys.executable, str(Path(__file__).resolve().parent / "phase_trgt_by_haplotag.py"),
           "--trgt-vcf", str(fixture_dir / "trgt.vcf"),
           "--spanning-bam", str(fixture_dir / "spanning.sam"),
           "--haplotag-list", str(fixture_dir / "haplotag.tsv"),
           "--allele-tag", "AL",
           "--out-alleles", str(fixture_dir / "alleles.tsv"),
           "--out-vcf", str(fixture_dir / "phased.vcf"),
           "--out-summary", str(out_summary), *extra]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(out_summary.read_text())


def allele_rows(fixture_dir: Path) -> dict[tuple[str, str], dict]:
    rows = {}
    text = (fixture_dir / "alleles.tsv").read_text().splitlines()
    header = text[0].split("\t")
    for line in text[1:]:
        row = dict(zip(header, line.split("\t")))
        rows[(row["locus_id"], row["hap"])] = row
    return rows


def vcf_records(fixture_dir: Path) -> dict[str, list[str]]:
    out = {}
    for line in (fixture_dir / "phased.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        trid = ptb.info_get(f[7], "TRID")
        out[trid] = f
    return out


def test_haplotag_list_whatshap_and_hiphase(tmp_path: Path):
    whatshap = tmp_path / "w.tsv"
    whatshap.write_text("#readname\thaplotype\tphaseset\tchromosome\n"
                        "a\tH1\t17020\tchr1\nb\tH2\t17020\tchr1\nc\tnone\t.\tchr1\n")
    haps, tally = ptb.read_haplotag_list(str(whatshap))
    assert haps == {"a": 1, "b": 2}
    assert tally["read_untagged"] == 1

    hiphase = tmp_path / "h.tsv"
    hiphase.write_text("source_block_index\tsample_name\tchrom\tphase_block_id\tread_name\thaplotag\n"
                       "0\tS1\tchr1\t17020\tr1\t2\n0\tS1\tchr1\t17020\tr2\t1\n")
    haps2, _ = ptb.read_haplotag_list(str(hiphase))
    assert haps2 == {"r1": 2, "r2": 1}


def test_reorder_sample_moves_every_per_slot_field():
    out = ptb.reorder_sample("GT:AL:ALLR:SD:MC:AP", "0/1:101,151:98-104,149-153:5,4:50,75:0.9,0.8",
                             [1, 0], phased=True)
    assert out == "1|0:151,101:149-153,98-104:4,5:75,50:0.8,0.9"


def test_reorder_sample_leaves_scalar_fields_alone():
    out = ptb.reorder_sample("GT:AL:DP", "0/1:101,151:30", [1, 0], phased=True)
    assert out == "1|0:151,101:30"


def test_assign_outcomes():
    perm, why, _ = ptb.assign([0, 0], {}, 2, 0.8)
    assert (perm, why) == ([0, 1], "homozygous")

    votes = {(0, 2): 5, (1, 1): 4}
    perm, why, support = ptb.assign([0, 1], votes, 2, 0.8)
    assert (perm, why) == ([1, 0], "reads")
    assert support[0] == (5, 1.0)

    perm, why, _ = ptb.assign([0, 1], {(0, 1): 5, (1, 1): 4}, 2, 0.8)
    assert (perm, why) == (None, "conflict")

    perm, why, _ = ptb.assign([0, 1], {(0, 1): 1, (1, 2): 1}, 2, 0.8)
    assert (perm, why) == (None, "insufficient_reads")

    # impure: 3 reads say hap1, 2 say hap2 -> purity 0.6 < 0.8
    perm, why, _ = ptb.assign([0, 1], {(0, 1): 3, (0, 2): 2, (1, 2): 4}, 2, 0.8)
    assert (perm, why) == (None, "insufficient_reads")


def test_end_to_end_resolution_tally(fixture_dir: Path):
    summary = run_script(fixture_dir)
    assert summary["status"] == "ok"
    assert summary["loci"] == 6
    assert summary["resolution"] == {
        "conflict": 1, "homozygous": 1, "insufficient_reads": 1, "reads": 3}
    # Five het loci are orientation-relevant; three of them resolved.
    assert summary["n_orientation_relevant"] == 5
    assert summary["frac_resolved_of_relevant"] == pytest.approx(0.6)
    assert summary["tag_census"]["AL"] == len(READS)
    assert summary["reads"]["read_joined"] == len(READS)


def test_swapped_locus_puts_long_allele_on_haplotype_one(fixture_dir: Path):
    run_script(fixture_dir)
    rows = allele_rows(fixture_dir)
    # reads place allele 1 (AL 151) on haplotype 1 at L_swap
    assert rows[("L_swap", "1")]["al_bp"] == "151"
    assert rows[("L_swap", "2")]["al_bp"] == "101"
    # and the unswapped locus keeps its order
    assert rows[("L_keep", "1")]["al_bp"] == "101"
    assert rows[("L_keep", "2")]["al_bp"] == "151"


def test_change_bp_and_units_use_reference_span(fixture_dir: Path):
    run_script(fixture_dir)
    rows = allele_rows(fixture_dir)
    # L_swap spans 1000-1100 inclusive = 101 bp; AL 151 is +50 bp = 25 AT units
    assert rows[("L_swap", "1")]["change_bp"] == "50"
    assert rows[("L_swap", "1")]["change_units"] == "25"
    assert rows[("L_swap", "2")]["change_bp"] == "0"


def test_unresolved_loci_are_absent_from_the_allele_table(fixture_dir: Path):
    run_script(fixture_dir)
    rows = allele_rows(fixture_dir)
    assert ("L_conflict", "1") not in rows
    assert ("L_thin", "1") not in rows
    # homozygous loci are emitted, flagged as such rather than as resolved phase
    assert rows[("L_hom", "1")]["phase_source"] == "homozygous"
    assert rows[("L_swap", "1")]["phase_source"] == "reads"


def test_vcf_phase_separator_and_allele_permutation(fixture_dir: Path):
    run_script(fixture_dir)
    recs = vcf_records(fixture_dir)
    swap = recs["L_swap"][9].split(":")
    assert swap[0] == "1|0", "resolved records are pipe-separated"
    assert swap[1] == "151,101", "AL must travel with GT"
    assert swap[3] == "75,50", "MC must travel with GT"
    assert ptb.info_get(recs["L_swap"][7], "TRGT_PHASE") == "reads"

    conflict = recs["L_conflict"][9].split(":")
    assert conflict[0] == "0/1", "unresolved records stay unphased"
    assert ptb.info_get(recs["L_conflict"][7], "TRGT_PHASE") == "conflict"
    assert ptb.info_get(recs["L_hom"][7], "TRGT_PHASE") == "homozygous"


def test_no_length_is_invented(fixture_dir: Path):
    """Every emitted pair is a permutation of the input pair, never a new value."""
    run_script(fixture_dir)
    recs = vcf_records(fixture_dir)
    original = {trid: al for _pos, _end, trid, _m, _gt, al, _sd, _mc in ROWS}
    for trid, f in recs.items():
        keys = f[8].split(":")
        al = f[9].split(":")[keys.index("AL")]
        assert sorted(al.split(",")) == sorted(original[trid].split(",")), trid


def test_missing_allele_tag_is_a_hard_error_with_the_census(fixture_dir: Path):
    cmd = [sys.executable, str(Path(__file__).resolve().parent / "phase_trgt_by_haplotag.py"),
           "--trgt-vcf", str(fixture_dir / "trgt.vcf"),
           "--spanning-bam", str(fixture_dir / "spanning.sam"),
           "--haplotag-list", str(fixture_dir / "haplotag.tsv"),
           "--allele-tag", "ZZ",
           "--out-summary", str(fixture_dir / "s.json")]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode != 0
    assert "--census-only" in proc.stderr


def test_census_only_needs_no_tag_or_haplotags(fixture_dir: Path):
    out = fixture_dir / "census.json"
    cmd = [sys.executable, str(Path(__file__).resolve().parent / "phase_trgt_by_haplotag.py"),
           "--trgt-vcf", str(fixture_dir / "trgt.vcf"),
           "--spanning-bam", str(fixture_dir / "spanning.sam"),
           "--census-only", "--out-summary", str(out)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    summary = json.loads(out.read_text())
    assert summary["status"] == "census_only"
    assert summary["tag_census"] == {
        "AL": len(READS), "FL": len(READS), "TR": len(READS)}


def test_samtools_failure_is_surfaced(tmp_path: Path):
    with pytest.raises(SystemExit) as exc:
        list(ptb.iter_bam(str(tmp_path / "missing.bam")))
    assert "samtools view failed" in str(exc.value)


def test_multiallelic_locus_is_oriented_by_slot_not_allele_index(fixture_dir: Path):
    """The bug this guards: the read tag is a genotype SLOT, not an ALT index.

    At L_multi the GT is 1/2, so the slots carry alleles 1 and 2 while the reads
    still carry AL:i:0 and AL:i:1. Keying votes by the allele index found nothing
    here and the locus fell through as insufficient_reads.
    """
    run_script(fixture_dir)
    rows = allele_rows(fixture_dir)
    # reads put slot 1 (AL 181) on haplotype 1 and slot 0 (AL 121) on haplotype 2
    assert rows[("L_multi", "1")]["al_bp"] == "181"
    assert rows[("L_multi", "2")]["al_bp"] == "121"
    assert rows[("L_multi", "1")]["phase_source"] == "reads"

    recs = vcf_records(fixture_dir)
    gt = recs["L_multi"][9].split(":")[0]
    assert gt == "2|1", "the ALT indices travel with the slots they label"


def test_slot_tag_out_of_range_is_counted_not_silently_used(fixture_dir: Path):
    """A tag holding a length rather than a slot must not be treated as a slot."""
    sam = fixture_dir / "spanning.sam"
    text = sam.read_text().replace("AL:i:0", "AL:i:101").replace("AL:i:1", "AL:i:151")
    sam.write_text(text)
    summary = run_script(fixture_dir)
    assert summary["reads"]["read_slot_out_of_range"] == len(READS)
    assert summary["reads"].get("read_joined", 0) == 0
    assert summary["resolution"].get("reads", 0) == 0


def test_locus_comes_from_the_tr_tag_by_default(fixture_dir: Path):
    summary = run_script(fixture_dir)
    assert summary["reads"]["read_locus_from_tag"] == len(READS)
    assert "read_ambiguous_locus" not in summary["reads"]


def test_overlap_fallback_still_works_when_the_tag_is_disabled(fixture_dir: Path):
    summary = run_script(fixture_dir, "--trid-tag", "")
    assert summary["reads"]["read_joined"] == len(READS)
    assert "read_locus_from_tag" not in summary["reads"]
    assert summary["resolution"]["reads"] == 3


def test_a_tagged_locus_absent_from_the_vcf_is_counted(fixture_dir: Path):
    sam = fixture_dir / "spanning.sam"
    sam.write_text(sam.read_text().replace("TR:Z:L_swap", "TR:Z:L_not_in_vcf"))
    summary = run_script(fixture_dir)
    assert summary["reads"]["read_locus_not_in_vcf"] == 4, "the four L_swap reads"
    assert "L_swap" not in allele_rows(fixture_dir).get(("L_swap", "1"), {})


def test_stopping_early_does_not_report_a_samtools_failure(fixture_dir: Path):
    """The census takes the first N records; SIGPIPE is not a failure.

    Closing the generator early made samtools exit non-zero, which the finally
    block raised as 'samtools view failed' with an empty stderr -- printed as
    'Exception ignored in: <generator object iter_bam>' after an otherwise
    successful run.
    """
    gen = ptb.iter_bam(str(fixture_dir / "spanning.sam"))
    first = next(gen)
    assert first[0].startswith("r_")
    gen.close()          # must not raise


def test_a_truncated_census_is_still_reported_honestly(fixture_dir: Path):
    summary = run_script(fixture_dir, "--census-records", "3")
    assert summary["tag_census_records"] == 3
    assert summary["status"] == "ok", "an early census must not abort the join"
    assert summary["reads"]["spanning_reads"] == len(READS), "the join still reads every record"


def test_census_defaults_to_the_whole_file(fixture_dir: Path):
    summary = run_script(fixture_dir)
    assert summary["tag_census_records"] == len(READS)
    assert summary["tag_census"]["TR"] == len(READS)

