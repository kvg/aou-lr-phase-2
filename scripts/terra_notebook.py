"""Terra Workbench notebook bootstrap — localize scripts/ from WORKSPACE_BUCKET."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def find_scripts_dir() -> Path | None:
    """Return local scripts/ if workspace_paths.py is present."""
    for d in (Path.cwd() / "scripts", Path.cwd().parent / "scripts"):
        if (d / "workspace_paths.py").is_file():
            return d.resolve()
    return None


def workspace_bucket() -> str:
    return os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")


def sync_scripts_from_bucket(dest: Path | None = None) -> Path:
    """Rsync the full scripts/ prefix from the workspace bucket."""
    bucket = workspace_bucket()
    if not bucket:
        raise FileNotFoundError(
            "scripts/ not found locally and WORKSPACE_BUCKET is unset. "
            "Upload repo scripts/ to gs://WORKSPACE/scripts/ or clone locally."
        )
    if dest is None:
        dest = Path.cwd() / "scripts"
    dest = dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    src = f"{bucket}/scripts/"
    print(f"gsutil -m rsync -r {src} {dest}/")
    subprocess.check_call(["gsutil", "-m", "rsync", "-r", src, str(dest) + "/"])
    return dest


def ensure_scripts_on_path(*required: str) -> Path:
    """Ensure scripts/ exists on sys.path and required files are present."""
    scripts_dir = find_scripts_dir()
    if scripts_dir is None:
        scripts_dir = sync_scripts_from_bucket()
    path = str(scripts_dir)
    if path not in sys.path:
        sys.path.insert(0, path)

    from workspace_paths import ensure_script_files

    if not required:
        required = ("workspace_paths.py",)
    return ensure_script_files(*required)


def init_notebook(*required: str) -> Path:
    """Notebook entry point: sync scripts/ if needed, verify CLIs, return scripts dir."""
    return ensure_scripts_on_path(*required)
