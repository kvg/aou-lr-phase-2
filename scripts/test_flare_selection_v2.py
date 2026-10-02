#!/usr/bin/env python3
"""Unit tests for recipe selection v2 (FLARE2 + whole-cohort LAI).

Covers: ancestry header / FLARE2-model AF mapping in the allele scorer,
bare-contig regions in the Mendelian scorer, within-population negative
controls, trio-bootstrap CIs and ``select_recipe``.
"""

from __future__ import annotations

import math
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_lai_exp import enrich_association_scores, select_recipe, trio_bootstrap_violation_ci  # noqa: E402
from flare_make_negative_control import derangement, within_population_mapping  # noqa: E402
from flare_score_allele_ancestry import (  # noqa: E402
    anc_weights_from_header,
    anc_weights_from_model,
    mixture_af,
    parse_ancestry_header,
    score_haplotype,
)
from flare_score_mendelian_lai import WHOLE_CONTIG_END, parse_region  # noqa: E402

FLARE2_MODEL = """\
anc0_afr\tanc1_eur
afr\teur\teas\tamr\tsas
10
0.5\t0.5
0.8\t0.2\t0\t0\t0
0.1\t0.9\t0\t0\t0
0.001\t0.001\t0.001\t0.001\t0.001
0.001\t0.001\t0.001\t0.001\t0.001
0.2\t0.3
"""


def test_ancestry_header_mapping():
    hdr = "##fileformat=VCFv4.2\n##ANCESTRY=<eas=0,amr=1,eur=2,afr=3,sas=4>\n#CHROM\tPOS\n"
    assert parse_ancestry_header(hdr) == {0: "eas", 1: "amr", 2: "eur", 3: "afr", 4: "sas"}
    w = anc_weights_from_header(hdr)
    assert w[3] == {"afr": 1.0}
    # Cluster names are not panels → must pass --model.
    try:
        anc_weights_from_header("##ANCESTRY=<anc0_afr=0,anc1_eur=1>\n")
    except SystemExit as exc:
        assert "--model" in str(exc)
    else:
        raise AssertionError("expected SystemExit for cluster ancestries without --model")


def test_flare2_model_mixture_af():
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "x.model"
        path.write_text(FLARE2_MODEL)
        w = anc_weights_from_model(path)
    assert w[0] == {"afr": 0.8, "eur": 0.2}
    afs = {"afr": 0.9, "eur": 0.1, "eas": 0.5, "amr": 0.5, "sas": 0.5}
    assert abs(mixture_af(w[0], afs) - (0.8 * 0.9 + 0.2 * 0.1)) < 1e-12
    ll, _, ok = score_haplotype(1, 0, afs, w)
    assert ok and abs(ll - math.log(0.74)) < 1e-12
    # Missing panel AF drops out and the rest renormalize.
    assert abs(mixture_af({"afr": 0.5, "eur": 0.5}, {"afr": 0.6, "eur": float("nan")}) - 0.6) < 1e-12
    # Default (no weights) keeps the legacy index map: 3 = afr.
    ll2, _, ok2 = score_haplotype(1, 3, {"afr": 0.9, "eur": 0.1})
    assert ok2 and abs(ll2 - math.log(0.9)) < 1e-12


def test_bare_contig_region():
    assert parse_region("chr20") == ("chr20", 1, WHOLE_CONTIG_END)
    assert parse_region("chr22:10-20") == ("chr22", 10, 20)


def test_negative_control_mapping():
    import random

    rng = random.Random(1)
    perm = derangement(list("abcdef"), rng)
    assert sorted(perm) == list("abcdef")
    assert all(x != y for x, y in zip("abcdef", perm))
    samples = [f"s{i}" for i in range(10)]
    pops = {s: ("AFR" if i < 6 else "EUR") for i, s in enumerate(samples)}
    pops["s9"] = "MID"  # singleton population keeps its own track
    m = within_population_mapping(samples, pops, seed=3)
    assert sorted(m.values()) == sorted(samples)
    assert all(pops[k] == pops[v] for k, v in m.items())
    assert m["s9"] == "s9"
    assert all(k != v for k, v in m.items() if pops[k] != "MID")
    assert within_population_mapping(samples, pops, seed=3) == m


def test_trio_bootstrap_ci():
    trios = [{"n_hard_violations": v, "n_informative": 100} for v in (1, 2, 3, 2, 1, 4)]
    ci = trio_bootstrap_violation_ci(trios)
    assert ci["ci_low"] <= 13 / 600 <= ci["ci_high"]
    assert trio_bootstrap_violation_ci(trios[:1]) == {}
    row = enrich_association_scores(
        {"experiment": "a"},
        mendel_json={"violations_per_informative_locus": 13 / 600, "trios": trios},
    )
    assert row["violations_ci_low"] == ci["ci_low"]


def _row(eid, viol, lo, hi, ll):
    return {
        "experiment": eid,
        "violations_per_informative_locus": viol,
        "violations_ci_low": lo,
        "violations_ci_high": hi,
        "mean_ll": ll,
        "excess_recomb_over_expected": 0.1,
    }


def test_select_recipe_winner_tie_and_invalid():
    ctl = {"a": "neg_a", "b": "neg_b"}
    rows = [
        _row("a", 0.010, 0.008, 0.012, -0.50),
        _row("b", 0.030, 0.025, 0.035, -0.51),
        _row("neg_a", 0.20, 0.18, 0.22, -0.70),
        _row("neg_b", 0.21, 0.19, 0.23, -0.71),
    ]
    d = select_recipe(rows, negative_controls=ctl)
    assert d["status"] == "winner" and d["winner"] == "a", d
    assert d["mendelian_separates_controls"] and d["mean_ll_separates_controls"]

    rows_tie = [*rows[:1], _row("b", 0.012, 0.009, 0.015, -0.51), *rows[2:]]
    d2 = select_recipe(rows_tie, negative_controls=ctl)
    assert d2["status"] == "tie_human_decision" and d2["tied"] == ["b"]
    assert d2["best_point_estimate"] == "a"

    # Control scores as well as the real recipe → metric cannot rank.
    rows_bad = [*rows[:2], _row("neg_a", 0.005, 0.004, 0.006, -0.70), rows[3]]
    d3 = select_recipe(rows_bad, negative_controls=ctl)
    assert d3["status"] == "metric_invalid_mendelian_does_not_beat_controls"
    assert d3["winner"] is None

    # Missing Mendelian fails closed; gate below threshold eliminates.
    rows_missing = [{**rows[0], "violations_per_informative_locus": None}, *rows[1:]]
    d4 = select_recipe(rows_missing, negative_controls=ctl)
    assert "mendelian_missing" in d4["eliminated"]["a"] and d4["winner"] == "b"
    d5 = select_recipe(rows, negative_controls=ctl, gates={"min_mean_ll": -0.505})
    assert d5["eliminated"]["b"] == ["mean_ll_below_gate"] and d5["winner"] == "a"
    assert Counter(d5["survivors"]) == Counter(["a"])


if __name__ == "__main__":
    test_ancestry_header_mapping()
    test_flare2_model_mixture_af()
    test_bare_contig_region()
    test_negative_control_mapping()
    test_trio_bootstrap_ci()
    test_select_recipe_winner_tie_and_invalid()
    print("ok")
