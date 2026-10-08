#!/usr/bin/env python3
"""Ancestry-alphabet handling in scripts/flare_score_mendelian_lai.py.

Regression for a silent failure: ``ANCESTRY`` hard-codes the five reference
panels, so a FLARE2 model with ``nanc`` 6 emitted code 5, every call with that
code was read as missing, and every locus touching it dropped out of
``n_informative_locus_calls``. The violation rate was then computed on a
non-random subset of loci and compared against ``nanc`` 5 recipes scored on all
of theirs.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCORER = ROOT / "scripts" / "flare_score_mendelian_lai.py"

pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

ANC_HEADER = ("##ANCESTRY=<anc0_afr=0,anc1_amr=1,anc2_eas=2,anc3_eas=3,"
              "anc4_eur=4,anc5_eur=5>")
LABELS = (
    "index\tname\tupstream_label\tdominant_panel\tdominant_weight\n"
    "0\tanc0_afr\tanc_1\tafr\t0.99\n"
    "1\tanc1_amr\tanc_4\tamr\t0.29\n"
    "2\tanc2_eas\tanc_3\teas\t0.99\n"
    "3\tanc3_eas\tanc_5\teas\t0.36\n"
    "4\tanc4_eur\tanc_2\teur\t0.70\n"
    "5\tanc5_eur\tanc_0\teur\t0.53\n"
)

# Four loci, all Mendelian-consistent (child hap1 from father, hap2 from
# mother). Two of them use ancestry code 5, which is outside the five-panel
# default. Per locus: (c1, c2, f1, f2, m1, m2).
LOCI = [
    (0, 1, 0, 2, 1, 3),
    (5, 1, 5, 2, 1, 3),   # touches code 5
    (2, 3, 2, 0, 3, 4),
    (4, 5, 4, 0, 5, 2),   # touches code 5
]


def _build(tmp: Path, with_header: bool) -> Path:
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr20,length=64444167>"]
    if with_header:
        lines.append(ANC_HEADER)
    lines += ['##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="a1">',
              '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="a2">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCH\tFA\tMO"]
    for i, (c1, c2, f1, f2, m1, m2) in enumerate(LOCI):
        lines.append("\t".join([
            "chr20", str((i + 1) * 1000), ".", "A", "G", ".", "PASS", ".", "AN1:AN2",
            f"{c1}:{c2}", f"{f1}:{f2}", f"{m1}:{m2}"]))
    plain = tmp / ("anc.hdr.vcf" if with_header else "anc.nohdr.vcf")
    plain.write_text("\n".join(lines) + "\n")
    gz = plain.with_suffix(".vcf.gz")
    subprocess.run(["bash", "-c", f"bgzip -cf {plain} > {gz} && bcftools index -tf {gz}"], check=True)
    (tmp / "trio.ped").write_text("ped,id,father,mother\nFAM,CH,FA,MO\n")
    (tmp / "chr20.map").write_text(
        "".join(f"chr20\trs{i}\t{(i + 1) * 1000 / 1e6 * 1.5:.6f}\t{(i + 1) * 1000}\n"
                for i in range(len(LOCI))))
    return gz


def _score(tmp: Path, gz: Path, *extra: str, out: str = "s.json"):
    r = subprocess.run(
        [sys.executable, str(SCORER), "--anc-vcf", str(gz), "--ped", str(tmp / "trio.ped"),
         "--map", str(tmp / "chr20.map"), "--region", "chr20", "--out", str(tmp / out), *extra],
        capture_output=True, text=True)
    res = json.loads((tmp / out).read_text()) if r.returncode == 0 else None
    return r, res


def test_ancestry_header_keeps_every_locus(tmp_path: Path):
    gz = _build(tmp_path, with_header=True)
    r, res = _score(tmp_path, gz)
    assert r.returncode == 0, r.stderr
    assert res["ancestry_alphabet"] == [0, 1, 2, 3, 4, 5]
    assert res["ancestry_alphabet_source"] == "##ANCESTRY header"
    # All four loci are informative; the two using code 5 are no longer dropped.
    assert res["n_informative_locus_calls"] == len(LOCI)
    assert res["n_hard_violations"] == 0
    assert res["grid_call_stats"].get("n_out_of_range_calls", 0) == 0


def test_out_of_range_codes_refuse_rather_than_shrink_denominator(tmp_path: Path):
    gz = _build(tmp_path, with_header=False)
    r, _ = _score(tmp_path, gz, out="bad.json")
    assert r.returncode != 0
    assert "outside [0, 1, 2, 3, 4]" in r.stderr
    assert "max code seen 5" in r.stderr


def test_allow_out_of_range_reproduces_the_old_silent_drop(tmp_path: Path):
    gz = _build(tmp_path, with_header=False)
    r, res = _score(tmp_path, gz, "--allow-out-of-range", out="old.json")
    assert r.returncode == 0, r.stderr
    # The two loci touching code 5 vanish from the denominator — the bug.
    assert res["n_informative_locus_calls"] == 2
    assert res["grid_call_stats"]["n_out_of_range_calls"] == 4
    assert res["grid_call_stats"]["max_code_seen"] == 5


def test_num_ancs_overrides_the_header(tmp_path: Path):
    gz = _build(tmp_path, with_header=True)
    r, res = _score(tmp_path, gz, "--num-ancs", "6", out="n6.json")
    assert r.returncode == 0, r.stderr
    assert res["ancestry_alphabet_source"] == "--num-ancs=6"
    assert res["n_informative_locus_calls"] == len(LOCI)


def test_projection_collapses_to_dominant_panels(tmp_path: Path):
    gz = _build(tmp_path, with_header=True)
    (tmp_path / "labels.tsv").write_text(LABELS)
    r, res = _score(tmp_path, gz, "--project-labels", str(tmp_path / "labels.tsv"), out="p.json")
    assert r.returncode == 0, r.stderr
    # eas, amr, eur, afr -> codes 0..3; the two eas and two eur clusters merge.
    assert res["projected_to"] == [0, 1, 2, 3]
    assert res["n_informative_locus_calls"] == len(LOCI)
    # Coarsening can only remove violations, never add them.
    assert res["n_hard_violations"] == 0


def test_projection_must_cover_every_ancestry(tmp_path: Path):
    gz = _build(tmp_path, with_header=True)
    short = "\n".join(LABELS.splitlines()[:5]) + "\n"   # drops codes 4 and 5
    (tmp_path / "short.tsv").write_text(short)
    r, _ = _score(tmp_path, gz, "--project-labels", str(tmp_path / "short.tsv"), out="sp.json")
    assert r.returncode != 0
    assert "does not cover ancestry codes" in r.stderr


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
