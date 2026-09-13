#!/usr/bin/env python3
"""Sync this repo onto an AoU Terra notebook VM (no laptop gsutil needed).

Designed to run **inside** a Terra Workbench notebook where ``WORKSPACE_BUCKET``
and Firecloud/FISS credentials already exist. Typical flow:

1. ``git clone`` / ``git pull`` ``kvg/aou-lr-phase-2`` (HTTPS + optional token)
2. ``gsutil -m rsync`` ``scripts/`` → ``$WORKSPACE_BUCKET/scripts/``
3. Copy ``notebooks/terra/*.ipynb`` onto the persistent notebook disk (+ optional bucket mirror)
4. Upload WDLs to ``$WORKSPACE_BUCKET/wdl/`` and optionally push new Firecloud method snapshots
5. Upsert data tables (default: ``flare_lai_exp`` from ``flare/configs/lai_exp.tsv``) via FISS

Example (Terra notebook)::

    from terra_sync_repo import sync_all
    sync_all(github_token=os.environ.get("GITHUB_TOKEN"))

CLI::

    python3 scripts/terra_sync_repo.py --ref main --upsert-tables flare_lai_exp
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

DEFAULT_REPO_URL = "https://github.com/kvg/aou-lr-phase-2.git"
DEFAULT_REF = "main"
DEFAULT_CLONE_DIR = "aou-lr-phase-2"

# Production WDLs we care about staging / optionally snapshotting.
DEFAULT_WDLS: tuple[str, ...] = (
    "flare/wdl/FlareByPopulation.wdl",
    "propagate_annotations/wdl/PropagateFlareAncestry.wdl",
    "propagate_annotations/wdl/AnnotateFlareGqDp.wdl",
    "propagate_annotations/wdl/PropagateAnnotations.wdl",
    "tractor_mix/wdl/TractorMixPilot.wdl",
    "tractor_mix/wdl/TractorMixGenome.wdl",
    "tractor_mix/wdl/SaigePilot.wdl",
    "felix/wdl/FelixPilot.wdl",
    "felix/wdl/FelixGenome.wdl",
    "sv_annotation/wdl/AnnotateSvCallset.wdl",
    "snv_stats/wdl/BcftoolsGlnexusStats.wdl",
    "methylation_stats/wdl/PbCpgSampleStats.wdl",
    "methylation_stats/wdl/PbCpgChromSites.wdl",
)

# entity_type -> repo-relative TSV (Terra flexible import; first column entity:…_id)
DEFAULT_TABLE_TSVS: dict[str, str] = {
    "flare_lai_exp": "flare/configs/lai_exp.tsv",
}

DEFAULT_NAMESPACE = "allofus-drc-wgs-LR-prodData"
DEFAULT_WORKSPACE = "AoU_DRC_LongReads_PhaseTwo_Storage"


@dataclass
class SyncReport:
    repo_dir: str = ""
    git_sha: str = ""
    git_ref: str = ""
    scripts_dest: str = ""
    notebooks_copied: list[str] = field(default_factory=list)
    wdls_staged: list[str] = field(default_factory=list)
    methods_updated: list[str] = field(default_factory=list)
    tables_upserted: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    dry_run: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def run(
    cmd: Sequence[str] | str,
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    if isinstance(cmd, str):
        print("+", cmd)
        return subprocess.run(
            cmd,
            shell=True,
            cwd=str(cwd) if cwd else None,
            check=check,
            text=True,
            env=env,
        )
    print("+", " ".join(map(str, cmd)))
    return subprocess.run(
        list(map(str, cmd)),
        cwd=str(cwd) if cwd else None,
        check=check,
        text=True,
        env=env,
    )


def workspace_bucket() -> str:
    return os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")


def terra_namespace_workspace() -> tuple[str, str]:
    ns = os.environ.get("TERRA_NAMESPACE", DEFAULT_NAMESPACE)
    ws = os.environ.get("TERRA_WORKSPACE", DEFAULT_WORKSPACE)
    # AoU sometimes exposes WORKSPACE_NAME as namespace/workspace
    raw = os.environ.get("WORKSPACE_NAME", "")
    if "/" in raw and (not os.environ.get("TERRA_NAMESPACE") or not os.environ.get("TERRA_WORKSPACE")):
        a, b = raw.split("/", 1)
        ns, ws = a, b
    return ns, ws


def _auth_repo_url(url: str, token: Optional[str]) -> str:
    if not token:
        return url
    # https://github.com/org/repo.git → https://x-access-token:TOKEN@github.com/org/repo.git
    m = re.match(r"https://([^/]+)/(.*)$", url.strip())
    if not m:
        return url
    host, rest = m.group(1), m.group(2)
    return f"https://x-access-token:{token}@{host}/{rest}"


def ensure_repo(
    *,
    dest: Path,
    repo_url: str = DEFAULT_REPO_URL,
    ref: str = DEFAULT_REF,
    github_token: Optional[str] = None,
    dry_run: bool = False,
) -> Path:
    """Clone or fetch+checkout ``ref`` into ``dest``. Returns repo root."""
    dest = dest.expanduser().resolve()
    token = github_token or os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    auth_url = _auth_repo_url(repo_url, token)

    if dry_run:
        print(f"[dry-run] ensure_repo dest={dest} ref={ref}")
        return dest

    if (dest / ".git").is_dir():
        run(["git", "remote", "set-url", "origin", auth_url], cwd=dest, check=False)
        run(["git", "fetch", "--tags", "--force", "origin"], cwd=dest)
        # Prefer origin/ref, then tag/sha
        rc = run(
            ["git", "rev-parse", "--verify", f"origin/{ref}"],
            cwd=dest,
            check=False,
        )
        if rc.returncode == 0:
            run(["git", "checkout", "-B", ref, f"origin/{ref}"], cwd=dest)
        else:
            run(["git", "checkout", "--detach", ref], cwd=dest)
            run(["git", "reset", "--hard", ref], cwd=dest)
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and any(dest.iterdir()):
            raise SystemExit(
                f"Clone destination {dest} exists and is not a git repo. "
                "Set AOU_LR_REPO_DIR to an empty path or existing clone."
            )
        run(
            [
                "git",
                "clone",
                "--branch",
                ref,
                "--single-branch",
                auth_url,
                str(dest),
            ],
            check=False,
        )
        if not (dest / ".git").is_dir():
            # branch might be a SHA; clone default then checkout
            run(["git", "clone", auth_url, str(dest)])
            run(["git", "checkout", "--detach", ref], cwd=dest)

    # Avoid leaking token into later `git remote -v` in shared logs when possible
    if token and repo_url.startswith("https://"):
        run(["git", "remote", "set-url", "origin", repo_url], cwd=dest, check=False)

    return dest


def git_describe(repo: Path) -> tuple[str, str]:
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=str(repo),
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    ref = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        cwd=str(repo),
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    return sha, ref


def stage_scripts(
    repo: Path,
    *,
    bucket: Optional[str] = None,
    dry_run: bool = False,
) -> str:
    """Rsync repo scripts/ to $WORKSPACE_BUCKET/scripts/ (Cromwell + notebook bootstrap)."""
    bucket = (bucket or workspace_bucket()).rstrip("/")
    if not bucket:
        raise SystemExit("WORKSPACE_BUCKET is unset; cannot stage scripts/")
    src = repo / "scripts"
    if not src.is_dir():
        raise SystemExit(f"missing {src}")
    dest = f"{bucket}/scripts"
    if dry_run:
        print(f"[dry-run] gsutil -m rsync -r {src}/ {dest}/")
        return dest
    run(["gsutil", "-m", "rsync", "-r", str(src) + "/", dest + "/"])
    return dest


def stage_notebooks(
    repo: Path,
    *,
    local_dest: Path,
    bucket: Optional[str] = None,
    mirror_to_bucket: bool = True,
    dry_run: bool = False,
) -> list[str]:
    """Copy notebooks/terra/*.ipynb onto the notebook disk (+ optional GCS mirror)."""
    src_dir = repo / "notebooks" / "terra"
    if not src_dir.is_dir():
        raise SystemExit(f"missing {src_dir}")
    local_dest = local_dest.expanduser().resolve()
    copied: list[str] = []
    if dry_run:
        print(f"[dry-run] copy {src_dir}/*.ipynb → {local_dest}/")
    else:
        local_dest.mkdir(parents=True, exist_ok=True)
        for nb in sorted(src_dir.glob("*.ipynb")):
            shutil.copy2(nb, local_dest / nb.name)
            copied.append(nb.name)
    bucket = (bucket or workspace_bucket()).rstrip("/")
    if mirror_to_bucket and bucket:
        dest = f"{bucket}/notebooks/terra"
        if dry_run:
            print(f"[dry-run] gsutil -m rsync -r {src_dir}/ {dest}/")
        else:
            run(["gsutil", "-m", "rsync", "-r", str(src_dir) + "/", dest + "/"])
    return copied


def stage_wdls(
    repo: Path,
    *,
    bucket: Optional[str] = None,
    wdl_paths: Sequence[str] = DEFAULT_WDLS,
    dry_run: bool = False,
) -> list[str]:
    """Upload selected WDLs to $WORKSPACE_BUCKET/wdl/<same relative path>."""
    bucket = (bucket or workspace_bucket()).rstrip("/")
    if not bucket:
        raise SystemExit("WORKSPACE_BUCKET is unset; cannot stage WDLs")
    staged: list[str] = []
    for rel in wdl_paths:
        src = repo / rel
        if not src.is_file():
            print(f"skip missing WDL: {rel}", file=sys.stderr)
            continue
        dest = f"{bucket}/wdl/{rel}"
        if dry_run:
            print(f"[dry-run] gsutil cp {src} {dest}")
        else:
            run(["gsutil", "cp", str(src), dest])
        staged.append(rel)
    return staged


def _firecloud_api() -> Any:
    from firecloud import api as fapi

    return fapi


def upsert_entities_tsv(
    tsv_path: Path,
    *,
    namespace: str,
    workspace: str,
    dry_run: bool = False,
) -> None:
    """Upsert a Terra data table from a flexible entity TSV (FISS)."""
    text = tsv_path.read_text()
    first = text.splitlines()[0] if text else ""
    if not first.startswith("entity:"):
        raise SystemExit(
            f"{tsv_path}: first column must be entity:<type>_id (got {first!r})"
        )
    if dry_run:
        n = max(0, len(text.splitlines()) - 1)
        print(f"[dry-run] upload_entities {namespace}/{workspace} rows≈{n} from {tsv_path}")
        return
    fapi = _firecloud_api()
    # Prefer flexible model when available (AoU / newer Firecloud).
    try:
        resp = fapi.upload_entities(namespace, workspace, text, model="flexible")
    except TypeError:
        resp = fapi.upload_entities(namespace, workspace, text)
    if getattr(resp, "status_code", 200) not in (200, 201):
        raise RuntimeError(
            f"upload_entities failed ({resp.status_code}): {getattr(resp, 'text', '')[:800]}"
        )
    print(f"upserted entities from {tsv_path.name} → {namespace}/{workspace}")


def update_repository_method(
    wdl_path: Path,
    *,
    method_namespace: str,
    method_name: str,
    synopsis: str,
    dry_run: bool = False,
) -> Optional[int]:
    """Push a new snapshot to the Firecloud Methods repo. Returns snapshot id if known."""
    wdl_text = wdl_path.read_text()
    if dry_run:
        print(
            f"[dry-run] update_repository_method {method_namespace}/{method_name} "
            f"from {wdl_path}"
        )
        return None
    fapi = _firecloud_api()
    # Historical FISS signatures vary; try common ones.
    errors: list[str] = []
    for kwargs in (
        {"namespace": method_namespace, "method": method_name, "synopsis": synopsis, "wdl": wdl_text},
        {"namespace": method_namespace, "name": method_name, "synopsis": synopsis, "payload": wdl_text},
    ):
        try:
            # positional form used by many firecloud builds:
            # update_repository_method(namespace, method, synopsis, wdl)
            resp = fapi.update_repository_method(
                method_namespace, method_name, synopsis, wdl_text
            )
            break
        except TypeError as exc:
            errors.append(str(exc))
            try:
                resp = fapi.update_repository_method(**kwargs)
                break
            except Exception as exc2:  # noqa: BLE001
                errors.append(str(exc2))
                resp = None
        except Exception as exc:  # noqa: BLE001
            errors.append(str(exc))
            resp = None
    else:
        resp = None
    if resp is None:
        raise RuntimeError(
            "update_repository_method failed; stage WDL to bucket and import in UI. "
            + "; ".join(errors[:3])
        )
    code = getattr(resp, "status_code", 200)
    if code not in (200, 201):
        raise RuntimeError(
            f"update_repository_method HTTP {code}: {getattr(resp, 'text', '')[:800]}"
        )
    snap = None
    try:
        body = resp.json()
        snap = body.get("snapshotId") or body.get("payloadObject", {}).get("snapshotId")
    except Exception:  # noqa: BLE001
        pass
    print(f"updated method {method_namespace}/{method_name} snapshot={snap}")
    return int(snap) if snap is not None else None


def upsert_configured_tables(
    repo: Path,
    *,
    tables: Optional[Sequence[str]] = None,
    namespace: Optional[str] = None,
    workspace: Optional[str] = None,
    dry_run: bool = False,
) -> list[str]:
    ns, ws = terra_namespace_workspace()
    namespace = namespace or ns
    workspace = workspace or ws
    wanted = list(tables) if tables is not None else list(DEFAULT_TABLE_TSVS)
    done: list[str] = []
    for name in wanted:
        rel = DEFAULT_TABLE_TSVS.get(name)
        if not rel:
            print(f"skip unknown table key {name!r}", file=sys.stderr)
            continue
        path = repo / rel
        if not path.is_file():
            print(f"skip missing table TSV: {rel}", file=sys.stderr)
            continue
        upsert_entities_tsv(path, namespace=namespace, workspace=workspace, dry_run=dry_run)
        done.append(name)
    return done


def sync_all(
    *,
    repo_url: str = DEFAULT_REPO_URL,
    ref: str = DEFAULT_REF,
    clone_dir: Optional[Path] = None,
    github_token: Optional[str] = None,
    notebook_dest: Optional[Path] = None,
    stage_scripts_flag: bool = True,
    stage_notebooks_flag: bool = True,
    stage_wdls_flag: bool = True,
    upsert_tables: Optional[Sequence[str]] = ("flare_lai_exp",),
    update_methods: Optional[Sequence[dict[str, str]]] = None,
    method_namespace: Optional[str] = None,
    mirror_notebooks_to_bucket: bool = True,
    dry_run: bool = False,
) -> SyncReport:
    """End-to-end Terra sync. Returns a JSON-serializable report."""
    report = SyncReport(dry_run=dry_run)
    clone_dir = Path(
        clone_dir
        or os.environ.get("AOU_LR_REPO_DIR")
        or (Path.cwd() / DEFAULT_CLONE_DIR)
    )
    repo = ensure_repo(
        dest=clone_dir,
        repo_url=repo_url,
        ref=ref,
        github_token=github_token,
        dry_run=dry_run,
    )
    report.repo_dir = str(repo)
    if not dry_run:
        sha, gref = git_describe(repo)
        report.git_sha = sha
        report.git_ref = gref
    else:
        report.git_ref = ref

    if stage_scripts_flag:
        report.scripts_dest = stage_scripts(repo, dry_run=dry_run)

    if stage_notebooks_flag:
        nb_dest = Path(
            notebook_dest
            or os.environ.get("AOU_LR_NOTEBOOK_DEST")
            or Path.cwd()
        )
        report.notebooks_copied = stage_notebooks(
            repo,
            local_dest=nb_dest,
            mirror_to_bucket=mirror_notebooks_to_bucket,
            dry_run=dry_run,
        )

    if stage_wdls_flag:
        report.wdls_staged = stage_wdls(repo, dry_run=dry_run)

    if upsert_tables:
        try:
            report.tables_upserted = upsert_configured_tables(
                repo, tables=upsert_tables, dry_run=dry_run
            )
        except Exception as exc:  # noqa: BLE001
            report.warnings.append(f"table upsert failed: {exc}")

    if update_methods:
        mns = method_namespace or os.environ.get("TERRA_METHOD_NAMESPACE") or ""
        for spec in update_methods:
            rel = spec["wdl"]
            name = spec.get("name") or Path(rel).stem
            synopsis = spec.get("synopsis") or f"AoU LR sync {report.git_sha[:12]}"
            ns = spec.get("namespace") or mns
            if not ns:
                report.warnings.append(
                    f"skip method update for {rel}: set TERRA_METHOD_NAMESPACE "
                    "or pass method_namespace="
                )
                continue
            try:
                snap = update_repository_method(
                    repo / rel,
                    method_namespace=ns,
                    method_name=name,
                    synopsis=synopsis,
                    dry_run=dry_run,
                )
                report.methods_updated.append(f"{ns}/{name}:{snap}")
            except Exception as exc:  # noqa: BLE001
                report.warnings.append(f"method update {ns}/{name}: {exc}")

    if report.wdls_staged and not report.methods_updated:
        report.warnings.append(
            "WDLs staged under $WORKSPACE_BUCKET/wdl/ — re-import in Terra "
            "Workflows UI (or pass update_methods=... / TERRA_METHOD_NAMESPACE) "
            "to push Firecloud method snapshots automatically."
        )
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-url", default=DEFAULT_REPO_URL)
    p.add_argument("--ref", default=os.environ.get("AOU_LR_REF", DEFAULT_REF))
    p.add_argument("--clone-dir", type=Path, default=None)
    p.add_argument("--notebook-dest", type=Path, default=None)
    p.add_argument("--github-token", default=None)
    p.add_argument("--skip-scripts", action="store_true")
    p.add_argument("--skip-notebooks", action="store_true")
    p.add_argument("--skip-wdls", action="store_true")
    p.add_argument(
        "--upsert-tables",
        nargs="*",
        default=["flare_lai_exp"],
        help="Entity types to upsert (keys in DEFAULT_TABLE_TSVS). Empty list skips.",
    )
    p.add_argument(
        "--update-method",
        action="append",
        default=[],
        help="WDL repo path to push as a Firecloud method snapshot (repeatable)",
    )
    p.add_argument("--method-namespace", default=None)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(argv)

    methods = [{"wdl": w} for w in args.update_method]
    report = sync_all(
        repo_url=args.repo_url,
        ref=args.ref,
        clone_dir=args.clone_dir,
        github_token=args.github_token,
        notebook_dest=args.notebook_dest,
        stage_scripts_flag=not args.skip_scripts,
        stage_notebooks_flag=not args.skip_notebooks,
        stage_wdls_flag=not args.skip_wdls,
        upsert_tables=args.upsert_tables,
        update_methods=methods or None,
        method_namespace=args.method_namespace,
        dry_run=args.dry_run,
    )
    print(json.dumps(report.to_dict(), indent=2))
    if args.report_json:
        args.report_json.write_text(json.dumps(report.to_dict(), indent=2) + "\n")
    return 0 if not report.warnings else 0  # warnings are non-fatal


if __name__ == "__main__":
    raise SystemExit(main())
