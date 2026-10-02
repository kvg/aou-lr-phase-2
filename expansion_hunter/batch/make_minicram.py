#!/usr/bin/env python3
"""Native Batch minicram worker (print-reads 0.1.2: python, no bash).

Downloads catalog, hg38, and the CRAI with google.cloud.storage. Leaves the
Nearline CRAM as gs:// so str-analysis IntervalReader range-fetches only the
catalog containers. Uploads minicram, crai, transfer stats, and this log.

With EH_TASK_ROOT set (the one-task minicram + genotype job), work happens in
EH_TASK_ROOT/t{BATCH_TASK_INDEX} on the host disk, so two tasks packed on one
VM never share files. The reference, catalog and minicram stay there for the
genotype container, and task.json / task.env carry this task's row to the
later runnables. If the minicram from an earlier attempt is already on GCS it
is downloaded instead of re-read from Nearline.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

from google.cloud import storage

VIP = "199.36.153.4"
HOSTS = (
    "restricted.googleapis.com",
    "storage.googleapis.com",
    "oauth2.googleapis.com",
    "www.googleapis.com",
    "accounts.google.com",
    "iamcredentials.googleapis.com",
)


def _start_monitor(work: Path):
    try:
        mon_cls = ResourceMonitor  # prepended by submit_batch._worker_source
    except NameError:
        try:
            from resource_monitor import ResourceMonitor as mon_cls
        except ImportError:
            return None
    mon = mon_cls(work)
    mon.start()
    return mon


def _dump_monitor(mon, work: Path, *, sample_id: str, stage: str, fused: bool = False) -> Path | None:
    if mon is None:
        return None
    stats = mon.stop()
    stats["sample_id"] = sample_id
    stats["stage"] = stage
    path = work / (f"{stage}.resources.tsv" if fused else "resource_stats.tsv")
    try:
        write_resource_tsv(path, stats)
    except NameError:
        from resource_monitor import write_resource_tsv as _write

        _write(path, stats)
    print(
        f"resource peak_rss_mib={stats['peak_rss_bytes'] / (1024**2):.1f} "
        f"cpu_cores_avg={stats['cpu_cores_avg']:.2f} wall_sec={stats['wall_sec']:.1f}",
        flush=True,
    )
    return path


def task_dir() -> Path:
    root = os.environ.get("EH_TASK_ROOT", "").strip()
    if not root:
        return Path("/work")
    return Path(root) / f"t{int(os.environ.get('BATCH_TASK_INDEX', '0'))}"


def env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise SystemExit(f"{name} is empty")
    return value


def apply_task_tsv() -> None:
    uri = os.environ.get("TASKS_TSV", "").strip()
    if not uri:
        return
    idx = int(os.environ.get("BATCH_TASK_INDEX", "0"))
    project = (
        os.environ.get("GCLOUD_PROJECT", "").strip()
        or os.environ.get("GOOGLE_CLOUD_PROJECT", "").strip()
    )
    if not project:
        raise SystemExit("GCLOUD_PROJECT is empty (needed to load TASKS_TSV)")
    work = task_dir()
    work.mkdir(parents=True, exist_ok=True)
    dest = work / "tasks.tsv"
    download(client(project), uri, dest, project)
    with dest.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if idx < 0 or idx >= len(rows):
        raise SystemExit(f"BATCH_TASK_INDEX {idx} out of range n={len(rows)}")
    row = {k: ("" if v is None else str(v)) for k, v in rows[idx].items() if k}
    for key, value in row.items():
        os.environ[key] = value
    print(f"TASKS_TSV row {idx} SAMPLE_ID={os.environ.get('SAMPLE_ID')}", flush=True)
    if os.environ.get("EH_TASK_ROOT", "").strip():
        write_task_files(work, row)


def write_task_files(work: Path, row: dict[str, str]) -> None:
    """task.json for the genotype container, task.env for host scripts."""
    (work / "task.json").write_text(json.dumps(row), encoding="utf-8")
    with (work / "task.env").open("w", encoding="utf-8") as fh:
        for key, value in row.items():
            fh.write(f"export {key}={json.dumps(value)}\n")


def split_gs(uri: str) -> tuple[str, str]:
    if not uri.startswith("gs://"):
        raise SystemExit(f"expected gs:// URI, got {uri}")
    bucket, _, blob = uri[5:].partition("/")
    if not bucket or not blob:
        raise SystemExit(f"bad gs:// URI: {uri}")
    return bucket, blob


def ensure_restricted_vip() -> None:
    try:
        existing = Path("/etc/hosts").read_text(encoding="utf-8")
    except OSError:
        return
    needed = [h for h in HOSTS if h not in existing]
    if not needed:
        return
    try:
        with Path("/etc/hosts").open("a", encoding="utf-8") as fh:
            for host in needed:
                fh.write(f"{VIP} {host}\n")
    except OSError:
        return


def client(project: str) -> storage.Client:
    os.environ.pop("GOOGLE_APPLICATION_CREDENTIALS", None)
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    os.environ["CLOUDSDK_CORE_PROJECT"] = project
    return storage.Client(project=project)


def download(gcs: storage.Client, uri: str, dest: Path, project: str) -> None:
    bucket_name, blob_name = split_gs(uri)
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"download {uri} -> {dest}", flush=True)
    gcs.bucket(bucket_name, user_project=project).blob(blob_name).download_to_filename(str(dest))
    if dest.stat().st_size <= 0:
        raise SystemExit(f"empty download: {uri}")


def blob_exists(gcs: storage.Client, uri: str, project: str) -> bool:
    bucket_name, blob_name = split_gs(uri)
    try:
        return gcs.bucket(bucket_name, user_project=project).blob(blob_name).exists()
    except Exception as exc:
        print(f"exists check failed for {uri}: {exc}", flush=True)
        return False


def upload(gcs: storage.Client, src: Path, uri: str, project: str) -> None:
    bucket_name, blob_name = split_gs(uri)
    print(f"upload {src} -> {uri} ({src.stat().st_size} bytes)", flush=True)
    gcs.bucket(bucket_name, user_project=project).blob(blob_name).upload_from_filename(str(src))


def run_minicram(work: Path, *, cram: str, crai: Path, catalog: Path, ref: Path, out_cram: Path, project: str) -> None:
    window = os.environ.get("WINDOW_SIZE", "1000")
    merge = os.environ.get("MERGE_REGIONS_DISTANCE", "1000")
    cmd = [
        sys.executable,
        "-u",
        "-m",
        "str_analysis.make_minicram_for_expansion_hunter",
        "-R",
        str(ref),
        "-c",
        str(catalog),
        "-i",
        str(crai),
        "-o",
        str(out_cram),
        "-w",
        window,
        "-d",
        merge,
        "--verbose",
        "--output-data-transfer-stats",
        "--gcloud-project",
        project,
        cram,
    ]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=work)


def main() -> int:
    ensure_restricted_vip()
    work = task_dir()
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    apply_task_tsv()

    sample_id = env("SAMPLE_ID")
    cram = env("CRAM")
    crai_uri = env("CRAI")
    ref_fa_uri = env("REF_FA")
    ref_fai_uri = env("REF_FAI")
    catalog_uri = env("CATALOG")
    project = env("GCLOUD_PROJECT")
    minicram_uri = env("MINICRAM")
    minicrai_uri = env("MINICRAI")
    stats_uri = env("TRANSFER_STATS")
    resource_uri = os.environ.get("RESOURCE_STATS", "").strip()
    log_uri = os.environ.get("WORKER_LOG_GCS", "").strip()
    max_retry = int(os.environ.get("MAX_RETRY", "3"))
    wait_time = int(os.environ.get("WAIT_TIME", "30"))
    monitor = _start_monitor(work)

    if not cram.startswith("gs://") or not crai_uri.startswith("gs://"):
        raise SystemExit(f"CRAM and CRAI must stay gs:// (not localized): {cram} {crai_uri}")

    fused = bool(os.environ.get("EH_TASK_ROOT", "").strip())
    reuse = os.environ.get("REUSE_MINICRAM", "1").strip() != "0"
    log_path = work / ("minicram.worker.log" if fused else "worker.log")
    log_fh = log_path.open("w", encoding="utf-8")

    class Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for stream in self.streams:
                stream.write(data)
                stream.flush()
            return len(data)

        def flush(self):
            for stream in self.streams:
                stream.flush()

    sys.stdout = Tee(sys.__stdout__, log_fh)
    sys.stderr = Tee(sys.__stderr__, log_fh)

    gcs = client(project)
    status = 0
    try:
        print(f"minicram start sample={sample_id} project={project}", flush=True)
        print(f"CRAM={cram}", flush=True)
        print(f"CRAI={crai_uri}", flush=True)

        catalog = work / "catalog.json"
        ref = work / "reference.fa"
        ref_fai = work / "reference.fa.fai"
        crai = work / f"{sample_id}.cram.crai"
        out_cram = work / f"{sample_id}.minicram.cram"

        download(gcs, catalog_uri, catalog, project)
        json.load(catalog.open())
        download(gcs, ref_fa_uri, ref, project)
        download(gcs, ref_fai_uri, ref_fai, project)

        if reuse and all(blob_exists(gcs, u, project) for u in (minicram_uri, minicrai_uri, stats_uri)):
            print("minicram already on GCS; downloading it instead of reading the CRAM", flush=True)
            download(gcs, minicram_uri, out_cram, project)
            download(gcs, minicrai_uri, Path(str(out_cram) + ".crai"), project)
            print("minicram done (reused)", flush=True)
            return 0

        download(gcs, crai_uri, crai, project)

        last_error = None
        for attempt in range(1, max_retry + 1):
            for leftover in work.glob("*.data_transfer_stats.tsv"):
                leftover.unlink()
            out_cram.unlink(missing_ok=True)
            Path(str(out_cram) + ".crai").unlink(missing_ok=True)
            try:
                run_minicram(
                    work,
                    cram=cram,
                    crai=crai,
                    catalog=catalog,
                    ref=ref,
                    out_cram=out_cram,
                    project=project,
                )
                last_error = None
                break
            except subprocess.CalledProcessError as exc:
                last_error = exc
                print(f"make_minicram failed attempt {attempt}/{max_retry}", flush=True)
                if attempt == max_retry:
                    raise
                time.sleep(wait_time)

        if last_error is not None:
            raise last_error
        out_crai = Path(str(out_cram) + ".crai")
        if not out_cram.is_file() or out_cram.stat().st_size <= 0:
            raise SystemExit("make_minicram wrote no CRAM")
        if not out_crai.is_file() or out_crai.stat().st_size <= 0:
            raise SystemExit("make_minicram wrote no CRAI")
        stats = sorted(work.glob("*.data_transfer_stats.tsv"))
        if not stats:
            raise SystemExit("make_minicram wrote no data_transfer_stats.tsv")

        upload(gcs, out_cram, minicram_uri, project)
        upload(gcs, out_crai, minicrai_uri, project)
        upload(gcs, stats[0], stats_uri, project)
        print("minicram done", flush=True)
    except Exception:
        status = 1
        traceback.print_exc()
    finally:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        stats_path = _dump_monitor(monitor, work, sample_id=sample_id, stage="minicram", fused=fused)
        log_fh.close()
        up = client(project)
        if stats_path is not None and resource_uri:
            try:
                upload(up, stats_path, resource_uri, project)
            except Exception as exc:
                print(f"resource stats upload failed: {exc}", file=sys.__stderr__)
        if log_uri:
            try:
                upload(up, log_path, log_uri, project)
            except Exception as exc:
                print(f"worker log upload failed: {exc}", file=sys.__stderr__)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
