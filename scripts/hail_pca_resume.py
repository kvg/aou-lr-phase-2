#!/usr/bin/env python3
"""Resume helpers for Hail PCA notebooks — detect completed stages on GCS/local."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any


def env_flag(name: str, *, default: bool = False) -> bool:
    value = os.environ.get(name, "")
    if not value:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def output_exists(path: str) -> bool:
    """Return True if a file or directory exists."""
    if not path:
        return False
    if path.startswith("gs://"):
        result = subprocess.run(
            ["gsutil", "-q", "stat", path.rstrip("/")],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    return Path(path).exists()


def hail_matrix_exists(path: str) -> bool:
    """Return True if a Hail MatrixTable path looks complete."""
    if not path:
        return False
    rows_marker = path.rstrip("/") + "/rows"
    if path.startswith("gs://"):
        result = subprocess.run(
            ["gsutil", "-q", "ls", rows_marker],
            capture_output=True,
            text=True,
        )
        return result.returncode == 0
    return (Path(path) / "rows").exists()


def read_json_uri(path: str) -> dict[str, Any] | None:
    if not output_exists(path):
        return None
    if path.startswith("gs://"):
        result = subprocess.run(
            ["gsutil", "cat", path],
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(result.stdout)
    return json.loads(Path(path).read_text())


def write_json_uri(path: str, payload: dict[str, Any]) -> None:
    text = json.dumps(payload, indent=2) + "\n"
    if path.startswith("gs://"):
        proc = subprocess.run(
            ["gsutil", "cp", "-", path],
            input=text,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"gsutil cp failed for {path}: {proc.stderr[:500]}")
        return
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text)


def params_match(saved: dict[str, Any] | None, current: dict[str, Any]) -> bool:
    if not saved:
        return False
    return saved == current


def describe_stage(
    label: str,
    path: str,
    *,
    force: bool,
    params_path: str | None = None,
    current_params: dict[str, Any] | None = None,
    matrix: bool = False,
) -> dict[str, Any]:
    exists = hail_matrix_exists(path) if matrix else output_exists(path)
    saved = read_json_uri(params_path) if params_path else None
    params_ok = params_match(saved, current_params or {}) if current_params else True
    if exists and not force and params_ok:
        action = "skip"
    elif exists and not force and not params_ok:
        action = "rerun (params changed)"
    else:
        action = "run"
    return {
        "stage": label,
        "path": path,
        "exists": exists,
        "force": force,
        "params_ok": params_ok,
        "action": action,
    }


def print_stage_plan(stages: list[dict[str, Any]]) -> None:
    print("Resume plan:")
    for stage in stages:
        flag = "SKIP" if stage["action"] == "skip" else "RUN"
        print(f"  [{flag}] {stage['stage']}: {stage['action']} ({stage['path']})")
