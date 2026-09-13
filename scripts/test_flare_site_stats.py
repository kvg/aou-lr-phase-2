#!/usr/bin/env python3
"""Unit tests for scripts/flare_site_stats.py (no bcftools required)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare_site_stats import parse_flare_retained_markers  # noqa: E402


SAMPLE_LOG = """
Start Time            :  Tue Jan 01 00:00:00 UTC 2026
Command line: java -jar flare.jar ...

Statistics
  reference samples :  2504
  target samples    :  812
  markers           :  57123

Wallclock Time      :  1 hour
"""


def test_parse_markers_line():
    assert parse_flare_retained_markers(SAMPLE_LOG) == 57123


def test_parse_markers_with_commas():
    text = "Statistics\n  markers           :  57,123\n"
    assert parse_flare_retained_markers(text) == 57123


def test_parse_missing():
    assert parse_flare_retained_markers("no stats here\n") is None


if __name__ == "__main__":
    test_parse_markers_line()
    test_parse_markers_with_commas()
    test_parse_missing()
    print("ok")
