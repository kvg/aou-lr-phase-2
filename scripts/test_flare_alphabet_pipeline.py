#!/usr/bin/env python3
"""The alphabet-comparability guard, exercised through the real Part 8 chain.

scorer CLI -> mendelian_lai_score.json -> enrich_association_scores -> select_recipe

The earlier unit test for the guard hand-built rows that already contained
``ancestry_alphabet`` / ``projected_to``. The real chain dropped those fields in
``enrich_association_scores``, so the guard never fired in the notebook. This test
starts from the scorer's own output so that cannot happen again.

Recipes (all on synthetic trios, one true demography):
  pin        five reference panels, no projection (the original-FLARE recipe)
  flare2_a/b FLARE2 clusters (the real nanc=6 model), projected with --project-model
Each has a label-noised negative control, as in the notebook.
"""

from __future__ import annotations

import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from flare_lai_exp import (  # noqa: E402
    _bcftools_any_record, backfill_mendelian_alphabet, enrich_association_scores,
    select_recipe,
)

SCORER = ROOT / "scripts" / "flare_score_mendelian_lai.py"
pytestmark = pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")

NANC6_MODEL = 'anc0_afr\tanc1_amr\tanc2_eas\tanc3_eas\tanc4_eur\tanc5_eur\neas\tamr\teur\tafr\tsas\n10\n0.16666667\t0.16666667\t0.16666667\t0.16666667\t0.16666667\t0.16666667\n0.0003265203\t4.7078331e-05\t0.00068987152\t0.99865249\t0.00028404098\n0.25814588\t0.29329949\t0.25031271\t0.080096667\t0.11814525\n0.99267025\t0.00046353228\t0.0023528348\t0.00031210448\t0.0042012835\n0.35752978\t0.0043622889\t0.32195037\t0.0039458728\t0.31221169\n0.0015301695\t0.0002097046\t0.70223458\t0.00047346685\t0.29555208\n0.023267487\t0.00042016028\t0.52837229\t0.30298554\t0.14495453\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\t9.3623174e-06\n430.76948\t959.38497\t340.2983\t709.65574\t303.34133\t926.54931\n'

PIN_HDR = "##ANCESTRY=<eas=0,amr=1,eur=2,afr=3,sas=4>"
F2_HDR = "##ANCESTRY=<anc0_afr=0,anc1_amr=1,anc2_eas=2,anc3_eas=3,anc4_eur=4,anc5_eur=5>"

# true panel -> FLARE2 cluster codes (the model's dominant panels: afr amr eas eas eur eur).
# FLARE2 has no SAS-dominant cluster, so SAS ancestry is expressed as EUR clusters.
PANEL_TO_CLUSTERS = {"eas": [2, 3], "amr": [1], "eur": [4, 5], "afr": [0], "sas": [4, 5]}
PIN_CODE = {"eas": 0, "amr": 1, "eur": 2, "afr": 3, "sas": 4}
PANELS = ["eas", "amr", "eur", "afr", "sas"]
WEIGHTS = [0.15, 0.15, 0.40, 0.25, 0.05]

N_TRIOS, N_LOCI = 24, 500


def _truth(rng: random.Random):
    """(child hap1, hap2, father hap1, hap2, mother hap1, hap2) true panels per locus, no error."""
    def mosaic():
        out, cur = [], None
        for i in range(N_LOCI):
            if cur is None or rng.random() < 0.02:
                cur = rng.choices(PANELS, WEIGHTS)[0]
            out.append(cur)
        return out
    f1, f2, m1, m2 = mosaic(), mosaic(), mosaic(), mosaic()
    return f1, f2, m1, m2


def _trio_labels(rng, to_code, eps, k_codes):
    f1, f2, m1, m2 = _truth(rng)
    c1, c2, flip = [], [], 0
    for i in range(N_LOCI):
        if rng.random() < 0.01:
            flip ^= 1
        c1.append((f1, f2)[flip][i]); c2.append((m1, m2)[flip][i])
    haps = [c1, c2, f1, f2, m1, m2]
    out = []
    for h in haps:
        row = []
        for p in h:
            code = to_code(p, rng)
            if rng.random() < eps:
                code = rng.choice([c for c in k_codes if c != code])
            row.append(code)
        out.append(row)
    return out


def _write(tmp: Path, name: str, to_code, eps, k_codes, header: str, seed: int) -> Path:
    rng = random.Random(seed)
    samples, cols = [], []
    for t in range(N_TRIOS):
        ch1, ch2, fa1, fa2, mo1, mo2 = _trio_labels(rng, to_code, eps, k_codes)
        samples += [f"c{t}", f"f{t}", f"m{t}"]
        cols += [(ch1, ch2), (fa1, fa2), (mo1, mo2)]
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr20,length=64444167>", header,
             '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="a1">',
             '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="a2">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for i in range(N_LOCI):
        cells = [f"{a[i]}:{b[i]}" for a, b in cols]
        lines.append("\t".join(["chr20", str((i + 1) * 1000), ".", "A", "G", ".", "PASS", ".", "AN1:AN2", *cells]))
    plain = tmp / f"{name}.vcf"
    plain.write_text("\n".join(lines) + "\n")
    gz = tmp / f"{name}.vcf.gz"
    subprocess.run(["bash", "-c", f"bgzip -cf {plain} > {gz} && bcftools index -tf {gz}"], check=True)
    return gz


def _shared_files(tmp: Path):
    (tmp / "trios.ped").write_text("ped,id,father,mother\n" + "".join(
        f"F{t},c{t},f{t},m{t}\n" for t in range(N_TRIOS)))
    (tmp / "chr20.map").write_text("".join(
        f"chr20\trs{i}\t{(i + 1) * 1000 / 1e6 * 1.5:.6f}\t{(i + 1) * 1000}\n" for i in range(N_LOCI)))
    (tmp / "nanc6.model").write_text(NANC6_MODEL)
    # SAS merged into EUR: pin projected to the four panels the FLARE2 clusters can express.
    (tmp / "pin_sas_eur.labels.tsv").write_text(
        "index\tname\tupstream_label\tdominant_panel\tdominant_weight\n"
        "0\teas\teas\teas\t1\n1\tamr\tamr\tamr\t1\n2\teur\teur\teur\t1\n3\tafr\tafr\tafr\t1\n4\tsas\tsas\teur\t1\n")


def _score(tmp: Path, name: str, gz: Path, *extra: str) -> dict:
    out = tmp / f"{name}.mendel.json"
    r = subprocess.run(
        [sys.executable, str(SCORER), "--anc-vcf", str(gz), "--ped", str(tmp / "trios.ped"),
         "--map", str(tmp / "chr20.map"), "--region", "chr20", "--out", str(out),
         "--experiment", name, *extra], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    return json.loads(out.read_text())


def _row(name: str, mendel: dict, mean_ll: float) -> dict:
    return enrich_association_scores({"experiment": name}, allele_json={"mean_ll": mean_ll},
                                     mendel_json=mendel)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("alpha")
    _shared_files(tmp)
    pin = lambda p, rng: PIN_CODE[p]  # noqa: E731
    f2 = lambda p, rng: rng.choice(PANEL_TO_CLUSTERS[p])  # noqa: E731
    vcfs = {
        "pin": _write(tmp, "pin", pin, 0.02, range(5), PIN_HDR, 1),
        "negctl_pin": _write(tmp, "negctl_pin", pin, 0.35, range(5), PIN_HDR, 2),
        "flare2_a": _write(tmp, "flare2_a", f2, 0.02, range(6), F2_HDR, 3),
        "negctl_flare2_a": _write(tmp, "negctl_flare2_a", f2, 0.35, range(6), F2_HDR, 4),
        "flare2_b": _write(tmp, "flare2_b", f2, 0.03, range(6), F2_HDR, 5),
        "negctl_flare2_b": _write(tmp, "negctl_flare2_b", f2, 0.35, range(6), F2_HDR, 6),
    }
    return tmp, vcfs


NEG = {"pin": "negctl_pin", "flare2_a": "negctl_flare2_a", "flare2_b": "negctl_flare2_b"}


def _rows(tmp, vcfs, pin_extra=(), strip_pin=False):
    model = ["--project-model", str(tmp / "nanc6.model")]
    rows = []
    for name, gz in vcfs.items():
        base = name.replace("negctl_", "")
        extra = model if base.startswith("flare2") else list(pin_extra)
        m = _score(tmp, name, gz, *extra)
        if strip_pin and base == "pin":
            for k in ("ancestry_alphabet", "projected_to", "n_ancestries"):
                m.pop(k, None)               # what a scorer predating those fields wrote
        rows.append(_row(name, m, -0.95 if name.startswith("negctl") else -0.50))
    return rows


def test_enrich_carries_the_alphabet_fields(world):
    tmp, vcfs = world
    row = _row("flare2_a", _score(tmp, "flare2_a_x", vcfs["flare2_a"], "--project-model",
                                   str(tmp / "nanc6.model")), -0.5)
    assert row["ancestry_alphabet"] == [0, 1, 2, 3, 4, 5]
    assert row["projected_to"] == [0, 1, 2, 3]


def test_pin_on_five_panels_vs_projected_flare2_is_refused(world):
    """The configuration the notebook ran: pin unprojected (5 labels), FLARE2 projected (4)."""
    tmp, vcfs = world
    d = select_recipe(_rows(tmp, vcfs), negative_controls=NEG)
    assert d["status"] == "metric_incomparable_mixed_ancestry_alphabets", d["status"]
    assert d["alphabets_comparable"] is False and d["winner"] is None
    assert d["ancestry_alphabets"]["pin"] == [0, 1, 2, 3, 4]
    assert d["ancestry_alphabets"]["flare2_a"] == [0, 1, 2, 3]


def test_pin_projected_to_the_same_four_panels_is_rankable(world):
    tmp, vcfs = world
    rows = _rows(tmp, vcfs, pin_extra=["--project-labels", str(tmp / "pin_sas_eur.labels.tsv")])
    d = select_recipe(rows, negative_controls=NEG)
    assert d["alphabets_comparable"] is True, d["ancestry_alphabets"]
    assert d["status"] in {"winner", "tie_human_decision"}, d["status"]
    assert {tuple(v) for v in d["ancestry_alphabets"].values()} == {(0, 1, 2, 3)}


def test_a_recipe_with_no_alphabet_info_is_not_assumed_to_match(world):
    """A pin scored by an older scorer has no alphabet fields; it must not pass as 'comparable'."""
    tmp, vcfs = world
    rows = _rows(tmp, vcfs, strip_pin=True)
    d = select_recipe(rows, negative_controls=NEG)
    assert d["status"] == "metric_incomparable_mixed_ancestry_alphabets", d["status"]
    assert d["alphabets_unknown"] and set(d["alphabets_unknown"]) <= {"pin", "negctl_pin"}


# ---------------------------------------------------------------------------
# Decision: keep the original-FLARE pin on five labels and rank it against the
# FLARE2 recipes projected to four, with an explicit override.
# ---------------------------------------------------------------------------

def _set_rates(rows, **rates):
    """Fix candidate rates (CI +/-1.5%) so the label-count rule is tested deterministically.

    Alphabets, controls and every other field stay as the real scorer wrote them."""
    out = []
    for r in rows:
        r = dict(r)
        if r["experiment"] in rates:
            v = rates[r["experiment"]]
            r.update(violations_per_informative_locus=v, violations_ci_low=v * 0.985,
                     violations_ci_high=v * 1.015)
        out.append(r)
    return out


def test_mixed_alphabets_still_refused_unless_the_override_is_explicit(world):
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs), pin=0.056, flare2_a=0.040, flare2_b=0.070)
    assert select_recipe(rows, negative_controls=NEG)["status"] == \
        "metric_incomparable_mixed_ancestry_alphabets"
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["mixed_alphabets_override"] is True and d["allow_mixed_alphabets"] is True
    assert d["alphabet_sizes"] == {"pin": 5, "flare2_a": 4, "flare2_b": 4}
    assert "caveat" in d


def test_fewer_label_winner_within_the_label_bias_is_a_tie_not_a_win(world):
    """flare2_a (4 labels) beats the pin (5) by 7.5%, less than one label's allowance."""
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs), pin=0.0430, flare2_a=0.0400, flare2_b=0.0600)
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["status"] == "tie_human_decision" and d["winner"] is None, d["status"]
    assert d["best_point_estimate"] == "flare2_a"
    assert d["tied"] == ["pin"] and d["tie_reason"] == "within_label_count_bias"
    cmp_pin = d["label_count_comparison_vs_best"]["pin"]
    assert cmp_pin["bias"] == "favours_best" and cmp_pin["within_bias"] is True
    assert cmp_pin["labels_best"] == 4 and cmp_pin["labels_competitor"] == 5
    # flare2_b shares flare2_a's alphabet, so the label count does not excuse the margin
    assert d["label_count_comparison_vs_best"]["flare2_b"]["bias"] == "none"
    assert "flare2_b" not in d["tied"]


def test_fewer_label_winner_clearing_the_label_bias_is_declared(world):
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs), pin=0.0560, flare2_a=0.0400, flare2_b=0.0700)  # 40% gap
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["status"] == "winner" and d["winner"] == "flare2_a", d["status"]
    assert d["label_count_comparison_vs_best"]["pin"]["within_bias"] is False


def test_a_pin_win_on_more_labels_is_conservative_and_stands(world):
    """The pin carries the label-count penalty, so beating FLARE2 by 5% (CIs disjoint) is a real win."""
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs), pin=0.0300, flare2_a=0.0315, flare2_b=0.0400)
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["status"] == "winner" and d["winner"] == "pin", d["status"]
    assert d["label_count_comparison_vs_best"]["flare2_a"]["bias"] == "conservative_for_best"


def test_unknown_alphabet_under_the_override_is_never_declared_a_winner(world):
    """If the pin's alphabet cannot be established the label-count effect is unbounded."""
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs, strip_pin=True), pin=0.0600, flare2_a=0.0300, flare2_b=0.0500)
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["status"] == "tie_human_decision" and "pin" in d["tied"], d["status"]
    assert d["label_count_comparison_vs_best"]["pin"]["bias"] == "unknown"


def test_the_override_does_not_bypass_the_negative_control_check(world):
    tmp, vcfs = world
    rows = _set_rates(_rows(tmp, vcfs), pin=0.0560, flare2_a=0.0400, flare2_b=0.0700,
                      negctl_flare2_a=0.0390)         # control no worse than its candidate
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["status"] == "metric_invalid_mendelian_does_not_beat_controls", d["status"]


def test_backfill_reads_the_pin_alphabet_from_the_vcf_header(world):
    """The pin's cached score predates the alphabet fields; the header recovers them."""
    tmp, vcfs = world
    m = _score(tmp, "pin_old", vcfs["pin"])
    for k in ("ancestry_alphabet", "ancestry_alphabet_source", "projected_to", "n_ancestries"):
        m.pop(k, None)
    fixed = backfill_mendelian_alphabet(m, str(vcfs["pin"]))
    assert fixed["ancestry_alphabet"] == [0, 1, 2, 3, 4] and fixed["alphabet_backfilled"] is True
    assert fixed["projected_to"] is None and fixed["violations_per_informative_locus"] == \
        m["violations_per_informative_locus"]
    again = backfill_mendelian_alphabet(fixed, str(vcfs["pin"]))
    assert again is fixed                            # idempotent: nothing left to fill
    # end to end: with the backfilled alphabet the override records 5 vs 4 labels
    rows = _rows(tmp, vcfs, strip_pin=True)
    for r in rows:
        if r["experiment"] in ("pin", "negctl_pin"):
            r.update(ancestry_alphabet=[0, 1, 2, 3, 4], n_ancestries=5)
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["alphabet_sizes"]["pin"] == 5 and d["alphabets_unknown"] == []


def test_backfill_refuses_when_the_old_scorer_would_have_dropped_loci(world):
    """FLARE2 VCFs declare codes 0-5; an old scorer silently dropped code 5, so its rate is void."""
    tmp, vcfs = world
    m = _score(tmp, "f2_old", vcfs["flare2_a"], "--project-model", str(tmp / "nanc6.model"))
    for k in ("ancestry_alphabet", "projected_to", "n_ancestries"):
        m.pop(k, None)
    with pytest.raises(ValueError, match=r"outside 0-4"):
        backfill_mendelian_alphabet(m, str(vcfs["flare2_a"]))


def test_backfill_refuses_without_an_ancestry_header(world):
    tmp, _ = world
    gz = _write(tmp, "nohdr", lambda p, rng: PIN_CODE[p], 0.02, range(5), "##source=nohdr", 9)
    with pytest.raises(ValueError, match=r"no ##ANCESTRY header"):
        backfill_mendelian_alphabet({"violations_per_informative_locus": 0.05}, str(gz))


def test_a_fewer_label_lead_is_flagged_provisional_and_a_more_label_lead_is_not(world):
    tmp, vcfs = world
    base = _rows(tmp, vcfs)
    clear = select_recipe(_set_rates(base, pin=0.0560, flare2_a=0.0400, flare2_b=0.0700),
                          negative_controls=NEG, allow_mixed_alphabets=True)
    assert clear["status"] == "winner" and clear["winner"] == "flare2_a"
    assert clear["provisional"] is True and clear["provisional_against"] == ["pin"]
    assert "PIN_SAS_TO" in clear["provisional_reason"]
    held = select_recipe(_set_rates(base, pin=0.0300, flare2_a=0.0315, flare2_b=0.0400),
                         negative_controls=NEG, allow_mixed_alphabets=True)
    assert held["winner"] == "pin" and held["provisional"] is False


def test_comparable_alphabets_are_never_provisional(world):
    tmp, vcfs = world
    rows = _rows(tmp, vcfs, pin_extra=["--project-labels", str(tmp / "pin_sas_eur.labels.tsv")])
    d = select_recipe(rows, negative_controls=NEG, allow_mixed_alphabets=True)
    assert d["mixed_alphabets_override"] is False and "provisional" not in d


# ---------------------------------------------------------------------------
# The real original-FLARE pin VCF has no ##ANCESTRY header (observed on the VM).
# ---------------------------------------------------------------------------

def _panel(tmp: Path) -> Path:
    path = tmp / "markers.tsv"
    path.write_text("chrom\tpos\n" + "".join(f"chr20\t{(i + 1) * 1000}\n" for i in range(N_LOCI)))
    return path


def _headerless(tmp, name, k_codes, seed=21):
    return _write(tmp, name, lambda p, rng: PIN_CODE[p], 0.02, k_codes, "##source=headerless", seed)


def test_headerless_pin_alphabet_is_established_from_the_data(world):
    tmp, _ = world
    gz = _headerless(tmp, "hl_pin", range(5))
    old = {"violations_per_informative_locus": 0.05}
    out = backfill_mendelian_alphabet(old, str(gz), panel=_panel(tmp))
    assert out["ancestry_alphabet"] == [0, 1, 2, 3, 4] and out["alphabet_backfilled"] is True
    assert "sampled" in out["ancestry_alphabet_source"] and "no ##ANCESTRY header" in out["ancestry_alphabet_source"]
    assert out["violations_per_informative_locus"] == 0.05      # the cached rate is carried unchanged


def test_headerless_with_a_code_above_four_is_refused(world):
    tmp, _ = world
    gz = _headerless(tmp, "hl_six", range(6))                   # FLARE2-like: codes 0-5 occur
    with pytest.raises(ValueError, match=r"above 4"):
        backfill_mendelian_alphabet({"x": 1}, str(gz), panel=_panel(tmp))


def test_headerless_that_never_shows_a_panel_code_is_not_confirmed(world):
    tmp, _ = world
    gz = _write(tmp, "hl_four", lambda p, rng: PIN_CODE[p] % 4, 0.02, range(4), "##source=h4", 22)
    with pytest.raises(ValueError, match=r"never occur"):
        backfill_mendelian_alphabet({"x": 1}, str(gz), panel=_panel(tmp))


def test_a_failed_bcftools_read_is_not_mistaken_for_a_clean_scan(tmp_path):
    with pytest.raises(RuntimeError, match="bcftools failed"):
        _bcftools_any_record(str(tmp_path / "missing.vcf.gz"), "chr20:1-100", "FMT/AN1>4")


def test_sampling_is_not_exhaustive_and_the_current_scorer_is(world):
    """One out-of-range call far from every sampled window slips past the sample, which is
    why the docstring calls it a sample; re-scoring with the current scorer does refuse it."""
    import gzip
    tmp, _ = world
    plain = tmp / "hl_stray.vcf"
    with gzip.open(_headerless(tmp, "hl_stray_src", range(5), seed=23), "rt") as fh:
        lines = fh.read().splitlines()
    first = next(i for i, l in enumerate(lines) if not l.startswith("#"))
    row = lines[first + 60].split("\t")                          # locus 60: outside every window
    row[9] = "5:5"
    lines[first + 60] = "\t".join(row)
    plain.write_text("\n".join(lines) + "\n")
    gz = tmp / "hl_stray.vcf.gz"
    subprocess.run(["bash", "-c", f"bgzip -cf {plain} > {gz} && bcftools index -tf {gz}"], check=True)
    out = backfill_mendelian_alphabet({"x": 1}, str(gz), panel=_panel(tmp), n_windows=5, window_bp=2000)
    assert out["ancestry_alphabet"] == [0, 1, 2, 3, 4]          # the sample misses the stray call
    r = subprocess.run(
        [sys.executable, str(SCORER), "--anc-vcf", str(gz), "--ped", str(tmp / "trios.ped"),
         "--map", str(tmp / "chr20.map"), "--region", "chr20", "--out", str(tmp / "stray.json")],
        capture_output=True, text=True)
    assert r.returncode != 0 and "outside" in (r.stderr + r.stdout)
