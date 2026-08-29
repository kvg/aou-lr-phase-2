"""Resolve Terra-vs-local roots for AoU-LR notebooks.

On Terra Workbench, CLIs live at ``$WORKSPACE_BUCKET/scripts/``. Notebooks call
``terra_notebook.init_notebook(...)`` to rsync or copy them locally before use.
Local git checkouts keep ``notebooks/terra/`` (and ``notebooks/rw/``) under
``notebooks/``, with ``scripts/`` as a sibling of ``notebooks/``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
WORKSPACE = SCRIPTS.parent


def data_root() -> Path:
    """Covariates, ``resources/``, and summaries.

    Local git checkout: ``tractor_mix/``. Terra: the workspace root, unless
    ``AOU_DATA_ROOT`` is set.
    """
    env = os.environ.get("AOU_DATA_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    local = WORKSPACE / "tractor_mix"
    if local.is_dir():
        return local
    return WORKSPACE


def sv_output_dir() -> Path:
    """AnnotateSvCallset downloads / figure outputs."""
    pkg = WORKSPACE / "sv_annotation"
    if pkg.is_dir():
        out = pkg / "outputs"
    else:
        out = WORKSPACE / "sv_outputs"
    out.mkdir(parents=True, exist_ok=True)
    return out


def add_scripts_to_path() -> Path:
    path = str(SCRIPTS)
    if path not in sys.path:
        sys.path.insert(0, path)
    return SCRIPTS


def ensure_script_files(*names: str) -> Path:
    """Verify named files under scripts/; fetch any missing ones from the bucket."""
    scripts_dir = SCRIPTS
    bucket = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")
    missing = [name for name in names if not (scripts_dir / name).is_file()]
    if not missing:
        return scripts_dir
    if not bucket:
        raise FileNotFoundError(
            f"Missing scripts: {', '.join(missing)}. "
            "Set WORKSPACE_BUCKET or copy scripts/ locally."
        )
    for name in missing:
        src = f"{bucket}/scripts/{name}"
        dest = scripts_dir / name
        print(f"gsutil cp {src} {dest}")
        subprocess.check_call(["gsutil", "cp", src, str(dest)])
    return scripts_dir
