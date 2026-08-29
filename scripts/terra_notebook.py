"""Terra Workbench notebook bootstrap — localize scripts/ from WORKSPACE_BUCKET."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def find_scripts_dir() -> Path | None:
    """Return local scripts/ if workspace_paths.py is present."""
    here = Path.cwd()
    for d in (here / "scripts", here.parent / "scripts", here.parent.parent / "scripts"):
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


def ensure_scripts_on_path(*required: str, refresh: bool | None = None) -> Path:
    """Ensure scripts/ exists on sys.path and required files are present."""
    if refresh is None:
        refresh = os.environ.get("TERRA_SYNC_SCRIPTS", os.environ.get("PCA_SYNC_SCRIPTS", "")).lower() in {
            "1",
            "true",
            "yes",
        }
    scripts_dir = find_scripts_dir()
    if scripts_dir is None or refresh:
        dest = scripts_dir if scripts_dir is not None else None
        scripts_dir = sync_scripts_from_bucket(dest)
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
