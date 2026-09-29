"""Thin gcloud storage wrappers (Jupyter / PET SA)."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=check, text=True, capture_output=True)


def exists(uri: str) -> bool:
    proc = _run(["gcloud", "storage", "ls", uri], check=False)
    return proc.returncode == 0


def cat(uri: str) -> str:
    proc = _run(["gcloud", "storage", "cat", uri])
    return proc.stdout


def upload_text(uri: str, text: str) -> None:
    proc = subprocess.run(
        ["gcloud", "storage", "cp", "-", uri],
        input=text,
        check=True,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr)


def upload_file(local: Path, uri: str) -> None:
    _run(["gcloud", "storage", "cp", str(local), uri])


def download_file(uri: str, local: Path) -> Path:
    local.parent.mkdir(parents=True, exist_ok=True)
    _run(["gcloud", "storage", "cp", uri, str(local)])
    return local


def ls(prefix: str) -> list[str]:
    proc = _run(["gcloud", "storage", "ls", prefix], check=False)
    if proc.returncode != 0:
        return []
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def ls_recursive(prefix: str) -> list[str]:
    proc = _run(["gcloud", "storage", "ls", "--recursive", prefix], check=False)
    if proc.returncode != 0:
        return []
    return [line.strip() for line in proc.stdout.splitlines() if line.strip().startswith("gs://")]


def object_size(uri: str) -> int | None:
    proc = _run(["gcloud", "storage", "ls", "--json", uri], check=False)
    if proc.returncode != 0 or not proc.stdout.strip():
        return None
    try:
        rows = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None
    if not rows:
        return None
    meta = rows[0].get("metadata") or rows[0]
    size = meta.get("size")
    return int(size) if size is not None else None
