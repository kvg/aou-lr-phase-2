#!/usr/bin/env python3
"""Tests for flare_build_indel_flanks command construction / BED span math."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "flare_build_indel_flanks.py"


def _load():
    spec = importlib.util.spec_from_file_location("flare_build_indel_flanks", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_indel_query_uses_type_filter_not_dash_v(monkeypatch, tmp_path: Path) -> None:
    mod = _load()
    calls: list[list[str]] = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            calls.append(cmd)
            self.stdout = iter(
                [
                    "chr20\t100\tA\tAT\n",  # insertion
                    "chr20\t200\tATGC\tA\n",  # deletion
                ]
            )
            self.stderr = None
            self.returncode = 0

        def communicate(self):
            return ("", "")

    monkeypatch.setattr(mod.subprocess, "Popen", FakeProc)
    rows = list(mod.iter_indel_intervals("x.vcf.gz", "chr20", flank=5))
    assert calls, "expected bcftools invoke"
    cmd = calls[0]
    assert cmd[:2] == ["bcftools", "query"]
    assert "-v" not in cmd
    assert "indels" not in cmd
    assert "-i" in cmd
    assert 'TYPE="indel"' in cmd
    assert "-r" in cmd and "chr20" in cmd
    # insertion at POS=100 REF=A ALT=AT, flank 5 → [94, 106)
    assert rows[0] == ("chr20", 94, 106)
    # deletion POS=200 REF=ATGC → end 203, +flank → [194, 208)
    assert rows[1] == ("chr20", 194, 208)


def test_bcftools_failure_surfaces_exit(monkeypatch) -> None:
    mod = _load()

    class FakeProc:
        def __init__(self, *a, **k):
            self.stdout = iter([])
            self.returncode = 255

        def communicate(self):
            return ("", "Error: failed to read indels\n")

    monkeypatch.setattr(mod.subprocess, "Popen", FakeProc)
    with pytest.raises(SystemExit) as ei:
        list(mod.iter_indel_intervals("x.vcf.gz", "", 5))
    assert "255" in str(ei.value)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
