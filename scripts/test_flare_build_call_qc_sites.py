#!/usr/bin/env python3
"""Unit tests for flare_build_call_qc_sites (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_build_call_qc_sites import rnc_has_i, sample_fails  # noqa: E402


def test_sample_fails_low_gq():
    assert sample_fails("5", "30", "..", min_gq=20, min_dp=10, drop_rnc_i=True) is True


def test_sample_fails_low_dp():
    assert sample_fails("40", "2", "..", min_gq=20, min_dp=10, drop_rnc_i=True) is True


def test_sample_fails_rnc_i():
    assert sample_fails("40", "30", "I.", min_gq=20, min_dp=10, drop_rnc_i=True) is True
    assert sample_fails("40", "30", "I.", min_gq=20, min_dp=10, drop_rnc_i=False) is False


def test_sample_passes():
    assert sample_fails("40", "30", "..", min_gq=20, min_dp=10, drop_rnc_i=True) is False


def test_sample_unqcable():
    assert sample_fails(".", ".", "..", min_gq=20, min_dp=10, drop_rnc_i=True) is None


def test_rnc_has_i():
    assert rnc_has_i("II")
    assert rnc_has_i("I.")
    assert not rnc_has_i("..")
    assert not rnc_has_i(".")


if __name__ == "__main__":
    test_sample_fails_low_gq()
    test_sample_fails_low_dp()
    test_sample_fails_rnc_i()
    test_sample_passes()
    test_sample_unqcable()
    test_rnc_has_i()
    print("ok")
