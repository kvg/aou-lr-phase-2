#!/usr/bin/env python3
"""Tests for terra_sync_repo helpers (no network / gsutil required)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from terra_sync_repo import (  # noqa: E402
    DEFAULT_TABLE_TSVS,
    DEFAULT_WDLS,
    SyncReport,
    _auth_repo_url,
)


def test_auth_repo_url():
    url = "https://github.com/kvg/aou-lr-phase-2.git"
    assert _auth_repo_url(url, None) == url
    authed = _auth_repo_url(url, "secret-token")
    assert authed.startswith("https://x-access-token:secret-token@github.com/")
    assert authed.endswith("/kvg/aou-lr-phase-2.git")


def test_defaults_point_at_real_paths():
    root = Path(__file__).resolve().parents[1]
    for rel in DEFAULT_WDLS:
        assert (root / rel).is_file(), rel
    for rel in DEFAULT_TABLE_TSVS.values():
        assert (root / rel).is_file(), rel
    tsv = (root / DEFAULT_TABLE_TSVS["flare_lai_exp"]).read_text().splitlines()[0]
    assert tsv.startswith("entity:flare_lai_exp_id")


def test_report_jsonable():
    r = SyncReport(repo_dir="/tmp/x", git_sha="abc", warnings=["w"])
    d = r.to_dict()
    assert d["repo_dir"] == "/tmp/x"
    assert d["warnings"] == ["w"]


if __name__ == "__main__":
    test_auth_repo_url()
    test_defaults_point_at_real_paths()
    test_report_jsonable()
    print("ok")
