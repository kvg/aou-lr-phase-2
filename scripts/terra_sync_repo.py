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
    notebooks_dest: str = ""
    notebooks_copied: list[str] = field(default_factory=list)
    wdls_staged: list[str] = field(default_factory=list)
    wdls_changed: list[str] = field(default_factory=list)
    wdls_unchanged: list[str] = field(default_factory=list)
    wdls_new: list[str] = field(default_factory=list)
    wdls_manual_import: list[str] = field(default_factory=list)
    methods_updated: list[str] = field(default_factory=list)
    configs_bumped: list[str] = field(default_factory=list)
    tables_upserted: list[str] = field(default_factory=list)
    submissions: list[str] = field(default_factory=list)
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
    """Clone or fetch+hard-reset ``ref`` into ``dest``. Returns repo root.

    The clone is a disposable mirror of GitHub. Local edits (e.g. Jupyter
    autosave inside ``CLONE_DIR``) are discarded on every sync.
    """
    dest = dest.expanduser().resolve()
    token = github_token or os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    auth_url = _auth_repo_url(repo_url, token)

    if dry_run:
        print(f"[dry-run] ensure_repo dest={dest} ref={ref}")
        return dest

    if (dest / ".git").is_dir():
        run(["git", "remote", "set-url", "origin", auth_url], cwd=dest, check=False)
        run(["git", "fetch", "--tags", "--force", "origin"], cwd=dest)
        rc = run(
            ["git", "rev-parse", "--verify", f"origin/{ref}"],
            cwd=dest,
            check=False,
        )
        target = f"origin/{ref}" if rc.returncode == 0 else ref
        # Discard dirty Jupyter autosaves / prior notebook copies in the mirror.
        run(["git", "reset", "--hard", target], cwd=dest)
        run(["git", "clean", "-fd"], cwd=dest)
        # Ensure branch name tracks ref when target is origin/ref
        if rc.returncode == 0:
            run(["git", "checkout", "-B", ref, target], cwd=dest, check=False)
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
            run(["git", "reset", "--hard", ref], cwd=dest)

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


def _copy_file(src: Path, dst: Path) -> None:
    """Copy file contents without requiring utime/xattr (Terra disks often forbid them)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copyfile(src, dst)
    except PermissionError:
        # Overwriting an open / root-owned notebook can fail; write via temp + replace.
        tmp = dst.with_suffix(dst.suffix + f".tmp.{os.getpid()}")
        try:
            shutil.copyfile(src, tmp)
            os.replace(tmp, dst)
        finally:
            if tmp.exists():
                tmp.unlink(missing_ok=True)


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
    skip_names: Sequence[str] = ("00_sync_repo.ipynb",),
) -> tuple[list[str], str]:
    """Copy notebooks/terra/*.ipynb to disk and ``$WORKSPACE_BUCKET/notebooks/``.

    Never writes into ``repo/notebooks/`` (that would dirty the disposable clone).
    Skips ``00_sync_repo.ipynb`` by default so a running sync notebook is not
    overwritten mid-execution.
    """
    src_dir = repo / "notebooks" / "terra"
    if not src_dir.is_dir():
        raise SystemExit(f"missing {src_dir}")
    local_dest = local_dest.expanduser().resolve()
    repo_resolved = repo.expanduser().resolve()
    skip = {Path(n).name for n in skip_names}
    copied: list[str] = []

    # If NOTEBOOK_DEST is inside the clone, stage beside the clone instead.
    try:
        local_dest.relative_to(repo_resolved)
        alt = repo_resolved.parent / "notebooks"
        print(
            f"NOTEBOOK_DEST {local_dest} is inside clone {repo_resolved}; "
            f"using {alt} instead to keep the git mirror clean",
            file=sys.stderr,
        )
        local_dest = alt
    except ValueError:
        pass

    if dry_run:
        print(f"[dry-run] copy {src_dir}/*.ipynb → {local_dest}/ (skip={sorted(skip)})")
    else:
        local_dest.mkdir(parents=True, exist_ok=True)
        for nb in sorted(src_dir.glob("*.ipynb")):
            if nb.name in skip:
                print(f"skip staging self/bootstrap notebook: {nb.name}", file=sys.stderr)
                continue
            dest = local_dest / nb.name
            try:
                _copy_file(nb, dest)
                copied.append(nb.name)
            except PermissionError as exc:
                print(
                    f"warning: could not copy {nb.name} → {dest} ({exc}); "
                    "continuing with GCS stage",
                    file=sys.stderr,
                )
    bucket = (bucket or workspace_bucket()).rstrip("/")
    notebooks_gcs = ""
    if mirror_to_bucket and bucket:
        # User-facing layout: $WORKSPACE_BUCKET/notebooks/*.ipynb
        notebooks_gcs = f"{bucket}/notebooks"
        if dry_run:
            excl = "|".join(re.escape(n) + "$" for n in sorted(skip)) or "a^"
            print(f"[dry-run] gsutil -m rsync -r -x '{excl}' {src_dir}/ {notebooks_gcs}/")
        else:
            # Exclude bootstrap notebook so a running sync does not overwrite itself
            # on the bucket either (leave any existing copy untouched).
            excl = "|".join(re.escape(n) + "$" for n in sorted(skip))
            cmd = ["gsutil", "-m", "rsync", "-r"]
            if excl:
                cmd.extend(["-x", excl])
            cmd.extend([str(src_dir) + "/", notebooks_gcs + "/"])
            run(cmd)
    return copied, notebooks_gcs


def _md5_local(path: Path) -> str:
    import hashlib

    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _md5_gcs(uri: str) -> Optional[str]:
    """Return hex MD5 of a GCS object, or None if missing / unavailable."""
    proc = subprocess.run(
        ["gsutil", "hash", "-m", "-h", uri],
        check=False,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        return None
    # gsutil hash -h prints "Hash (md5):\t<base64 or hex depending on version>"
    # Prefer parsing md5 line; fall back to gsutil ls -L Hash (md5)
    for line in (proc.stdout or "").splitlines():
        if "md5" in line.lower():
            parts = line.split(":", 1)
            if len(parts) == 2:
                val = parts[1].strip()
                # Some gsutil builds print base64; decode to hex when needed
                if re.fullmatch(r"[0-9a-fA-F]{32}", val):
                    return val.lower()
                try:
                    import base64

                    raw = base64.b64decode(val)
                    if len(raw) == 16:
                        return raw.hex()
                except Exception:  # noqa: BLE001
                    continue
    proc2 = subprocess.run(
        ["gsutil", "ls", "-L", uri],
        check=False,
        text=True,
        capture_output=True,
    )
    if proc2.returncode != 0:
        return None
    for line in (proc2.stdout or "").splitlines():
        if "Hash (md5)" in line or "md5=" in line.lower():
            # e.g. "        Hash (md5):           abc=..."
            m = re.search(r"md5\)?:\s*(\S+)", line, re.I)
            if not m:
                m = re.search(r"md5=([^\s,]+)", line, re.I)
            if not m:
                continue
            val = m.group(1).strip().rstrip(",")
            if re.fullmatch(r"[0-9a-fA-F]{32}", val):
                return val.lower()
            try:
                import base64

                raw = base64.b64decode(val)
                if len(raw) == 16:
                    return raw.hex()
            except Exception:  # noqa: BLE001
                continue
    return None


def stage_wdls(
    repo: Path,
    *,
    bucket: Optional[str] = None,
    wdl_paths: Sequence[str] = DEFAULT_WDLS,
    dry_run: bool = False,
) -> dict[str, list[str]]:
    """Upload selected WDLs; return staged/changed/new/unchanged/manual_import lists.

    ``wdls_manual_import`` is the GCS URIs (or repo-relative paths) that differ
    from the previous bucket copy — these need a Terra Workflows UI refresh when
    Methods-repo Create is unavailable.
    """
    bucket = (bucket or workspace_bucket()).rstrip("/")
    if not bucket:
        raise SystemExit("WORKSPACE_BUCKET is unset; cannot stage WDLs")
    staged: list[str] = []
    changed: list[str] = []
    new: list[str] = []
    unchanged: list[str] = []
    manual: list[str] = []
    for rel in wdl_paths:
        src = repo / rel
        if not src.is_file():
            print(f"skip missing WDL: {rel}", file=sys.stderr)
            continue
        dest = f"{bucket}/wdl/{rel}"
        local_md5 = _md5_local(src)
        remote_md5 = None if dry_run else _md5_gcs(dest)
        status: str
        if remote_md5 is None:
            status = "new"
        elif remote_md5 == local_md5:
            status = "unchanged"
        else:
            status = "changed"

        if dry_run:
            print(f"[dry-run] gsutil cp {src} {dest}  ({status})")
        else:
            if status != "unchanged":
                run(["gsutil", "cp", str(src), dest])
            else:
                print(f"unchanged (skip upload): {rel}", file=sys.stderr)

        staged.append(rel)
        if status == "new":
            new.append(rel)
            manual.append(dest)
        elif status == "changed":
            changed.append(rel)
            manual.append(dest)
        else:
            unchanged.append(rel)
    return {
        "staged": staged,
        "changed": changed,
        "new": new,
        "unchanged": unchanged,
        "manual_import": manual,
    }


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
    if dry_run:
        print(
            f"[dry-run] update_repository_method {method_namespace}/{method_name} "
            f"from {wdl_path}"
        )
        return None
    fapi = _firecloud_api()
    wdl_text = wdl_path.read_text()
    errors: list[str] = []
    resp = None
    # Prefer path form (FISS docs say wdl is a file); fall back to string body.
    for call in (
        lambda: fapi.update_repository_method(
            method_namespace, method_name, synopsis, str(wdl_path)
        ),
        lambda: fapi.update_repository_method(
            method_namespace, method_name, synopsis, wdl_text
        ),
    ):
        try:
            resp = call()
            break
        except Exception as exc:  # noqa: BLE001
            errors.append(str(exc))
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


def get_workspace_config(
    *,
    namespace: str,
    workspace: str,
    config_namespace: str,
    config_name: str,
) -> dict[str, Any]:
    fapi = _firecloud_api()
    resp = fapi.get_workspace_config(namespace, workspace, config_namespace, config_name)
    if getattr(resp, "status_code", 200) != 200:
        raise RuntimeError(
            f"get_workspace_config failed ({resp.status_code}): "
            f"{getattr(resp, 'text', '')[:800]}"
        )
    return resp.json()


def bump_workspace_config_snapshot(
    *,
    namespace: str,
    workspace: str,
    config_namespace: str,
    config_name: str,
    method_namespace: str,
    method_name: str,
    snapshot_id: int,
    inputs_overlay: Optional[dict[str, Any]] = None,
    dry_run: bool = False,
) -> str:
    """Point an existing workspace method config at a new method snapshot."""
    label = f"{config_namespace}/{config_name}@{snapshot_id}"
    if dry_run:
        print(f"[dry-run] bump workspace config {label}")
        return label
    cfg = get_workspace_config(
        namespace=namespace,
        workspace=workspace,
        config_namespace=config_namespace,
        config_name=config_name,
    )
    mrm = dict(cfg.get("methodRepoMethod") or {})
    mrm["methodNamespace"] = method_namespace
    mrm["methodName"] = method_name
    mrm["methodVersion"] = int(snapshot_id)
    # Keep URI in sync when present (Agora-style).
    mrm["methodUri"] = (
        f"agora://{method_namespace}/{method_name}/{int(snapshot_id)}"
    )
    cfg["methodRepoMethod"] = mrm
    if inputs_overlay:
        inputs = dict(cfg.get("inputs") or {})
        for k, v in inputs_overlay.items():
            inputs[k] = v if isinstance(v, str) else json.dumps(v)
        cfg["inputs"] = inputs
    fapi = _firecloud_api()
    resp = fapi.overwrite_workspace_config(
        namespace, workspace, config_namespace, config_name, cfg
    )
    if getattr(resp, "status_code", 200) not in (200, 201):
        raise RuntimeError(
            f"overwrite_workspace_config failed ({resp.status_code}): "
            f"{getattr(resp, 'text', '')[:800]}"
        )
    print(f"bumped workspace config → {label}")
    return label


def script_inputs_for_bucket(bucket: str) -> dict[str, str]:
    """FlareByPopulation script File inputs → current workspace bucket."""
    b = bucket.rstrip("/")
    return {
        "FlareByPopulation.split_script": f"{b}/scripts/flare_split_samples.py",
        "FlareByPopulation.summarize_script": f"{b}/scripts/flare_summarize_models.py",
        "FlareByPopulation.flare_model_script": f"{b}/scripts/flare_model.py",
        "FlareByPopulation.flare_site_stats_script": f"{b}/scripts/flare_site_stats.py",
        "FlareByPopulation.indel_flank_script": f"{b}/scripts/flare_build_indel_flanks.py",
    }


def list_entity_names(
    *,
    namespace: str,
    workspace: str,
    entity_type: str,
) -> list[dict[str, Any]]:
    fapi = _firecloud_api()
    resp = fapi.get_entities(namespace, workspace, entity_type)
    if getattr(resp, "status_code", 200) != 200:
        raise RuntimeError(
            f"get_entities failed ({resp.status_code}): {getattr(resp, 'text', '')[:800]}"
        )
    ents = resp.json()
    if not isinstance(ents, list):
        raise RuntimeError(f"unexpected entities payload: {type(ents)}")
    return ents


def incomplete_flare_lai_exp_ids(
    entities: Sequence[dict[str, Any]],
) -> list[str]:
    """Rows missing MergePopulationFlare outputs (anc_vcf + models_tsv)."""
    out: list[str] = []
    for ent in entities:
        name = str(ent.get("name") or "")
        attrs = ent.get("attributes") or {}
        anc = attrs.get("anc_vcf")
        models = attrs.get("models_tsv")
        anc_s = "" if anc is None else str(anc)
        models_s = "" if models is None else str(models)
        if not (anc_s.startswith("gs://") and models_s.startswith("gs://")):
            out.append(name)
    return out


def create_submission(
    *,
    namespace: str,
    workspace: str,
    config_namespace: str,
    config_name: str,
    entity_type: str,
    entity_name: str,
    use_callcache: bool = True,
    dry_run: bool = False,
) -> str:
    if dry_run:
        msg = (
            f"[dry-run] create_submission {config_namespace}/{config_name} "
            f"on {entity_type}/{entity_name}"
        )
        print(msg)
        return msg
    fapi = _firecloud_api()
    resp = fapi.create_submission(
        namespace,
        workspace,
        config_namespace,
        config_name,
        entity=entity_name,
        etype=entity_type,
        use_callcache=use_callcache,
    )
    code = getattr(resp, "status_code", 200)
    if code not in (200, 201):
        raise RuntimeError(
            f"create_submission failed ({code}) for {entity_name}: "
            f"{getattr(resp, 'text', '')[:800]}"
        )
    sid = None
    try:
        sid = resp.json().get("submissionId")
    except Exception:  # noqa: BLE001
        pass
    label = f"{entity_name}:{sid or 'ok'}"
    print("submitted", label)
    return label


def submit_flare_lai_exp(
    *,
    namespace: Optional[str] = None,
    workspace: Optional[str] = None,
    config_namespace: Optional[str] = None,
    config_name: str = "FlareByPopulation",
    entity_ids: Optional[Sequence[str]] = None,
    only_incomplete: bool = True,
    use_callcache: bool = True,
    dry_run: bool = False,
) -> list[str]:
    """Submit FlareByPopulation for selected / incomplete flare_lai_exp rows."""
    ns, ws = terra_namespace_workspace()
    namespace = namespace or ns
    workspace = workspace or ws
    config_namespace = config_namespace or namespace
    ents = list_entity_names(
        namespace=namespace, workspace=workspace, entity_type="flare_lai_exp"
    )
    if entity_ids is not None:
        wanted = set(entity_ids)
        names = [e["name"] for e in ents if e.get("name") in wanted]
        missing = wanted - set(names)
        if missing:
            raise SystemExit(f"unknown flare_lai_exp ids: {sorted(missing)}")
    elif only_incomplete:
        names = incomplete_flare_lai_exp_ids(ents)
    else:
        names = [str(e.get("name")) for e in ents if e.get("name")]
    if not names:
        print("no flare_lai_exp rows to submit")
        return []
    print(f"submitting {len(names)} flare_lai_exp row(s): {names}")
    return [
        create_submission(
            namespace=namespace,
            workspace=workspace,
            config_namespace=config_namespace,
            config_name=config_name,
            entity_type="flare_lai_exp",
            entity_name=name,
            use_callcache=use_callcache,
            dry_run=dry_run,
        )
        for name in names
    ]


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
    bump_configs: Optional[Sequence[dict[str, Any]]] = None,
    submit_flare: bool = False,
    submit_entity_ids: Optional[Sequence[str]] = None,
    submit_only_incomplete: bool = True,
    submit_config_name: str = "FlareByPopulation",
    submit_config_namespace: Optional[str] = None,
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

    ns, ws = terra_namespace_workspace()
    bucket = workspace_bucket()

    if stage_scripts_flag:
        report.scripts_dest = stage_scripts(repo, dry_run=dry_run)

    if stage_notebooks_flag:
        nb_dest = Path(
            notebook_dest
            or os.environ.get("AOU_LR_NOTEBOOK_DEST")
            or Path.cwd()
        )
        copied, gcs = stage_notebooks(
            repo,
            local_dest=nb_dest,
            mirror_to_bucket=mirror_notebooks_to_bucket,
            dry_run=dry_run,
        )
        report.notebooks_copied = copied
        report.notebooks_dest = gcs

    if stage_wdls_flag:
        wdl_info = stage_wdls(repo, dry_run=dry_run)
        report.wdls_staged = wdl_info["staged"]
        report.wdls_changed = wdl_info["changed"]
        report.wdls_new = wdl_info["new"]
        report.wdls_unchanged = wdl_info["unchanged"]
        report.wdls_manual_import = wdl_info["manual_import"]

    if upsert_tables:
        try:
            report.tables_upserted = upsert_configured_tables(
                repo, tables=upsert_tables, dry_run=dry_run
            )
        except Exception as exc:  # noqa: BLE001
            report.warnings.append(f"table upsert failed: {exc}")

    snap_by_method: dict[str, int] = {}
    if update_methods:
        mns = method_namespace or os.environ.get("TERRA_METHOD_NAMESPACE") or ns
        for spec in update_methods:
            rel = spec["wdl"]
            name = spec.get("name") or Path(rel).stem
            synopsis = spec.get("synopsis") or f"AoU LR sync {report.git_sha[:12] or ref}"
            method_ns = spec.get("namespace") or mns
            if not method_ns:
                report.warnings.append(
                    f"skip method update for {rel}: set TERRA_METHOD_NAMESPACE "
                    "or pass method_namespace="
                )
                continue
            try:
                snap = update_repository_method(
                    repo / rel,
                    method_namespace=method_ns,
                    method_name=name,
                    synopsis=synopsis,
                    dry_run=dry_run,
                )
                report.methods_updated.append(f"{method_ns}/{name}:{snap}")
                if snap is not None:
                    snap_by_method[f"{method_ns}/{name}"] = int(snap)
            except Exception as exc:  # noqa: BLE001
                msg = str(exc)
                if "403" in msg or "Authorization exception" in msg:
                    report.warnings.append(
                        f"method update {method_ns}/{name}: 403 — no Create on Methods "
                        f"namespace {method_ns!r}. WDLs are still at "
                        f"$WORKSPACE_BUCKET/wdl/{rel}. In Terra: Workflows → "
                        f"{name} → ⋯ → Version / Import → use that GCS path "
                        f"(or set METHOD_NAMESPACE to an Agora namespace you own)."
                    )
                else:
                    report.warnings.append(f"method update {method_ns}/{name}: {exc}")

    if bump_configs:
        for spec in bump_configs:
            method_ns = spec.get("method_namespace") or method_namespace or ns
            method_name = spec.get("method_name") or spec.get("name") or "FlareByPopulation"
            key = f"{method_ns}/{method_name}"
            snap = spec.get("snapshot_id") or snap_by_method.get(key)
            if snap is None:
                report.warnings.append(
                    f"skip config bump {spec}: no snapshot_id (method update may have failed)"
                )
                continue
            overlay = spec.get("inputs_overlay")
            if overlay is None and bucket and method_name == "FlareByPopulation":
                overlay = script_inputs_for_bucket(bucket)
            try:
                label = bump_workspace_config_snapshot(
                    namespace=ns,
                    workspace=ws,
                    config_namespace=spec.get("config_namespace") or ns,
                    config_name=spec.get("config_name") or method_name,
                    method_namespace=method_ns,
                    method_name=method_name,
                    snapshot_id=int(snap),
                    inputs_overlay=overlay,
                    dry_run=dry_run,
                )
                report.configs_bumped.append(label)
            except Exception as exc:  # noqa: BLE001
                report.warnings.append(f"config bump failed: {exc}")

    if report.wdls_manual_import and not report.methods_updated:
        lines = "\n".join(f"  - {u}" for u in report.wdls_manual_import)
        report.warnings.append(
            "Manual Terra Workflows import needed for updated/new WDLs "
            "(Methods-repo Create usually 403 on AoU):\n" + lines
        )
    elif report.wdls_staged and not report.methods_updated and not report.wdls_manual_import:
        report.warnings.append(
            "All staged WDLs match the bucket copies — no Workflows UI re-import needed."
        )

    if submit_flare:
        try:
            report.submissions = submit_flare_lai_exp(
                namespace=ns,
                workspace=ws,
                config_namespace=submit_config_namespace or ns,
                config_name=submit_config_name,
                entity_ids=submit_entity_ids,
                only_incomplete=submit_only_incomplete,
                dry_run=dry_run,
            )
        except Exception as exc:  # noqa: BLE001
            report.warnings.append(f"workflow submit failed: {exc}")

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
    p.add_argument(
        "--bump-config",
        action="append",
        default=[],
        help="Workspace config name to point at the new method snapshot (repeatable)",
    )
    p.add_argument(
        "--submit-flare",
        action="store_true",
        help="Submit FlareByPopulation on incomplete (or listed) flare_lai_exp rows",
    )
    p.add_argument(
        "--submit-entity",
        action="append",
        default=[],
        help="flare_lai_exp id to submit (repeatable; default=all incomplete)",
    )
    p.add_argument("--submit-all", action="store_true", help="Submit all rows, not only incomplete")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(argv)

    methods = [{"wdl": w} for w in args.update_method]
    bumps = [{"config_name": c, "method_name": c} for c in args.bump_config]
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
        bump_configs=bumps or None,
        submit_flare=args.submit_flare,
        submit_entity_ids=args.submit_entity or None,
        submit_only_incomplete=not args.submit_all,
        dry_run=args.dry_run,
    )
    print(json.dumps(report.to_dict(), indent=2))
    if args.report_json:
        args.report_json.write_text(json.dumps(report.to_dict(), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
