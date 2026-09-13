#!/usr/bin/env python3
"""Unit tests for scripts/flare_model.py."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_model import (  # noqa: E402
    build_input_model,
    match_pop_model,
    parse_model,
    pop_of_model_path,
    rewrite_t,
    set_props,
    write_model,
)

FIXTURE_DIR = Path(__file__).resolve().parent / "testdata" / "flare_models"
EM_MODEL = FIXTURE_DIR / "em_all_pops.AFR.model"
PIN_MODEL = FIXTURE_DIR / "pin_all_pops_gen_by_pop.AFR.model"


def test_roundtrip_byte_identical():
    for path in (EM_MODEL, PIN_MODEL):
        assert path.is_file(), f"missing fixture {path}"
        model = parse_model(path)
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / path.name
            write_model(model, out, preserve_raw=True)
            assert out.read_bytes() == path.read_bytes(), path.name


def test_parse_em_blocks():
    m = parse_model(EM_MODEL)
    assert m.ancestries == ["eas", "amr", "eur", "afr", "sas"]
    assert abs(m.t_gen - 93.7237) < 1e-6
    assert abs(sum(m.props) - 1.0) < 1e-4
    assert abs(m.props[3] - 0.741161) < 1e-6  # afr
    assert len(m.panel_weights) == 5 and len(m.panel_weights[0]) == 5
    assert len(m.miscopy) == 5 and len(m.miscopy[0]) == 5
    assert len(m.ibd_rates) == 5
    assert abs(m.mean_miscopy()[0] - sum(m.miscopy[0]) / 5) < 1e-12


def test_parse_pin_uniform_props():
    m = parse_model(PIN_MODEL)
    assert abs(m.t_gen - 8.0) < 1e-9
    assert all(abs(x - 0.2) < 1e-9 for x in m.props)


def test_props_sum_raises():
    m = parse_model(PIN_MODEL)
    m.props = [0.5, 0.5, 0.5, 0.0, 0.0]
    try:
        m.validate()
    except ValueError as exc:
        assert "sum" in str(exc).lower()
    else:
        raise AssertionError("expected ValueError")


def test_set_props_and_rewrite_t():
    m = parse_model(EM_MODEL)
    m2 = rewrite_t(m, 8.0)
    assert abs(m2.t_gen - 8.0) < 1e-9
    m3 = set_props(m2, {"afr": 0.8, "eur": 0.2})
    assert abs(m3.props[m3.ancestries.index("afr")] - 0.8) < 1e-9
    assert abs(sum(m3.props) - 1.0) < 1e-9
    try:
        set_props(m, {"afr": 0.1}, min_sum=0.5)
    except ValueError as exc:
        assert "typo" in str(exc) or "0.5" in str(exc)
    else:
        raise AssertionError("expected typo failure")
    try:
        set_props(m, {"MID": 1.0})
    except ValueError as exc:
        assert "MID" in str(exc)
    else:
        raise AssertionError("expected unknown ancestry failure")


def test_pop_of_and_match():
    assert pop_of_model_path(EM_MODEL) == "AFR"
    assert pop_of_model_path("x.AMR.model") == "AMR"
    try:
        match_pop_model("AFR", [EM_MODEL, PIN_MODEL], require=True)
    except ValueError as exc:
        assert "multiple" in str(exc)
    else:
        raise AssertionError("expected multiple-match failure")
    got = match_pop_model("AFR", [EM_MODEL], require=True)
    assert got == EM_MODEL
    try:
        match_pop_model("EUR", [EM_MODEL], require=True)
    except ValueError as exc:
        assert "EUR" in str(exc)
    else:
        raise AssertionError("expected missing pop failure")


def test_build_input_model_override():
    model, sources = build_input_model(
        pop="AFR",
        gen=10.0,
        gen_by_pop="AFR:8,AMR:12",
        pop_model_path=EM_MODEL,
        template_model_path=None,
        props_by_pop="AFR:AFR=0.80,EUR=0.20",
        inherit_props=True,
        inherit_panel_weights=False,
        inherit_mu=False,
    )
    assert model is not None
    assert abs(model.t_gen - 8.0) < 1e-9
    assert sources["t_gen"] == "override"
    assert sources["props"] == "override"
    assert sources["panel_weights"] == "default"
    assert sources["miscopy"] == "default"
    assert abs(model.props[model.ancestries.index("afr")] - 0.8) < 1e-9


if __name__ == "__main__":
    test_roundtrip_byte_identical()
    test_parse_em_blocks()
    test_parse_pin_uniform_props()
    test_props_sum_raises()
    test_set_props_and_rewrite_t()
    test_pop_of_and_match()
    test_build_input_model_override()
    print("ok")
