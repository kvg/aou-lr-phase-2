"""Terra Workbench notebook bootstrap — localize scripts/ from WORKSPACE_BUCKET.

Prefer refreshing the bucket from GitHub via ``notebooks/terra/00_sync_repo.ipynb``
(``scripts/terra_sync_repo.py``) when working on AoU Research Program Terra
workspaces (no laptop ``gsutil`` to the bucket). This module then rsyncs
``$WORKSPACE_BUCKET/scripts/`` onto the notebook disk for imports.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
from pathlib import Path


def find_scripts_dir() -> Path | None:
    """Return local scripts/ if workspace_paths.py is present."""
    here = Path.cwd()
    for d in (here / "scripts", here.parent / "scripts", here.parent.parent / "scripts"):
        if (d / "workspace_paths.py").is_file() or (d / "terra_notebook.py").is_file():
            return d.resolve()
    return None


def workspace_bucket() -> str:
    return os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")


def env_flag(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def drop_cached_script_modules(scripts_dir: Path) -> None:
    """Drop in-memory copies so a kernel re-run sees files just pulled from GCS."""
    importlib.invalidate_caches()
    prefix = str(scripts_dir.resolve())
    for name, mod in list(sys.modules.items()):
        path = getattr(mod, "__file__", None)
        if not path:
            continue
        try:
            resolved = str(Path(path).resolve())
        except OSError:
            continue
        if resolved == prefix or resolved.startswith(prefix + os.sep):
            sys.modules.pop(name, None)


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
    drop_cached_script_modules(dest)
    return dest


def ensure_scripts_on_path(*required: str, refresh: bool | None = None) -> Path:
    """Ensure scripts/ exists on sys.path and required files are present.

    On Terra, a persistent ``edit/scripts/`` copy is often stale. When
    ``WORKSPACE_BUCKET`` is set, rsync from the bucket unless
    ``TERRA_SYNC_SCRIPTS=false``.
    """
    bucket = workspace_bucket()
    already = env_flag("TERRA_SCRIPTS_LOCALIZED", default=False)
    if refresh is None:
        if already:
            refresh = False
        else:
            refresh = env_flag("TERRA_SYNC_SCRIPTS", default=bool(bucket)) or env_flag(
                "PCA_SYNC_SCRIPTS", default=False
            )
    scripts_dir = find_scripts_dir()
    if bucket and refresh:
        dest = scripts_dir if scripts_dir is not None else Path.cwd() / "scripts"
        scripts_dir = sync_scripts_from_bucket(dest)
        os.environ["TERRA_SCRIPTS_LOCALIZED"] = "true"
    elif scripts_dir is None:
        scripts_dir = sync_scripts_from_bucket()
        os.environ["TERRA_SCRIPTS_LOCALIZED"] = "true"
    path = str(scripts_dir)
    if path not in sys.path:
        sys.path.insert(0, path)

    from workspace_paths import ensure_script_files

    if not required:
        required = ("workspace_paths.py",)
    try:
        return ensure_script_files(*required, refresh=bool(bucket and refresh))
    except TypeError:
        return ensure_script_files(*required)


def init_notebook(*required: str) -> Path:
    """Notebook entry point: sync scripts/ if needed, verify CLIs, return scripts dir."""
    return ensure_scripts_on_path(*required)
