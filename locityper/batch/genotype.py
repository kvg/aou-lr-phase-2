#!/usr/bin/env python3
"""Native Batch genotype worker (Locityper image: locityper + python3, no GCS).

Host runnable already put minicram, fasta, jellyfish counts, BED, and vcf_db
on /work. This process runs preproc + per-locus genotype + summary CSV, then
the host uploads tar.gz / CSV.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import tarfile
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise SystemExit(f"{name} is empty")
    return value


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


def _dump_monitor(mon, work: Path, *, sample_id: str, stage: str) -> Path | None:
    if mon is None:
        return None
    stats = mon.stop()
    stats["sample_id"] = sample_id
    stats["stage"] = stage
    path = work / "resource_stats.tsv"
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


def apply_task_file() -> None:
    path = Path("/work/task.json")
    if not path.is_file():
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    for key, value in data.items():
        if key:
            os.environ[key] = "" if value is None else str(value)


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


def bed_loci(bed: Path) -> list[str]:
    names: list[str] = []
    for raw in bed.read_text(encoding="utf-8").splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 4 or not parts[3].strip():
            raise SystemExit(f"BED line missing column 4 (locus name): {raw}")
        names.append(parts[3].strip())
    if not names:
        raise SystemExit(f"no loci in {bed}")
    return names


def genotype_locus(locus: str, work: Path) -> None:
    out = work / "out_dir"
    (out / "loci" / locus).mkdir(parents=True, exist_ok=True)
    cmd = [
        "locityper",
        "genotype",
        "-a",
        str(work / "subset.cram"),
        "-r",
        str(work / "reference.fa"),
        "-d",
        str(work / "vcf_db"),
        "-p",
        str(work / "locityper_preproc"),
        "--subset-loci",
        locus,
        "-o",
        str(out),
    ]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=work, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def write_summary(work: Path, sample_id: str, dest: Path) -> int:
    loci_root = work / "out_dir" / "loci"
    rows = ["sample\tlocus\tgenotype\tquality\ttotal_reads\tunexpl_reads\tweight_dist\twarnings"]
    if loci_root.is_dir():
        for entry in sorted(loci_root.iterdir(), key=lambda p: p.name):
            if not entry.is_dir():
                continue
            json_path = entry / "res.json.gz"
            line = f"{sample_id}\t{entry.name}\t"
            if not json_path.is_file():
                rows.append(line + "*")
                continue
            with gzip.open(json_path, "rt") as fh:
                res = json.load(fh)
            if "genotype" not in res:
                rows.append(line + "*")
                continue
            gt = res["genotype"]
            qual = math.floor(10 * float(res["quality"])) * 0.1
            total_reads = res.get("total_reads", "")
            unexpl_reads = res.get("unexpl_reads", "")
            weight_dist = res.get("weight_dist")
            weight_s = "" if weight_dist is None else f"{float(weight_dist):.5f}"
            warnings = ";".join(res.get("warnings", [])) or "*"
            rows.append(
                line + f"{gt}\t{qual:.1f}\t{total_reads}\t{unexpl_reads}\t{weight_s}\t{warnings}"
            )
    dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return len(rows) - 1


def main() -> int:
    work = Path("/work")
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    log_path = work / "worker.log"
    log_fh = log_path.open("w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, log_fh)
    sys.stderr = Tee(sys.__stderr__, log_fh)
    status = 0
    monitor = None
    sample_id = "unknown"
    try:
        apply_task_file()
        sample_id = env("SAMPLE_ID")
        monitor = _start_monitor(work)
        n_cpu = int(os.environ.get("LOCITYPER_N_CPU", "2"))
        technology = os.environ.get("TECHNOLOGY", "illumina").strip() or "illumina"

        reads = work / "reads.cram"
        reads_idx = work / "reads.cram.crai"
        ref = work / "reference.fa"
        bed = work / "loci.bed"
        counts = work / "counts.jf"
        db_tar = work / "vcf_db.tar.gz"
        for path in (reads, reads_idx, ref, work / "reference.fa.fai", bed, counts, db_tar):
            if not path.is_file() or path.stat().st_size <= 0:
                raise SystemExit(f"missing input: {path}")

        (work / "subset.cram").unlink(missing_ok=True)
        (work / "subset.cram.crai").unlink(missing_ok=True)
        (work / "subset.cram").symlink_to(reads)
        (work / "subset.cram.crai").symlink_to(reads_idx)

        if shutil.which("locityper") is None:
            raise SystemExit("locityper is not on PATH")

        print(f"genotype start sample={sample_id} n_cpu={n_cpu} tech={technology}", flush=True)
        preproc = [
            "locityper",
            "preproc",
            "-a",
            str(work / "subset.cram"),
            "-r",
            str(ref),
            "-j",
            str(counts),
            "-@",
            str(n_cpu),
            "--technology",
            technology,
            "-o",
            str(work / "locityper_preproc"),
        ]
        print("+", " ".join(preproc), flush=True)
        subprocess.run(preproc, check=True, cwd=work)

        with tarfile.open(db_tar, "r:gz") as tar:
            try:
                tar.extractall(work, filter="data")
            except TypeError:
                tar.extractall(work)
        if not (work / "vcf_db").is_dir():
            raise SystemExit("vcf_db.tar.gz did not contain vcf_db/")

        loci = bed_loci(bed)
        print(f"loci: {len(loci)}", flush=True)
        workers = max(1, min(n_cpu, len(loci)))
        errors: list[str] = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(genotype_locus, locus, work): locus for locus in loci}
            for fut in as_completed(futs):
                locus = futs[fut]
                try:
                    fut.result()
                except Exception as exc:
                    errors.append(f"{locus}: {exc}")
        if errors:
            raise SystemExit("genotype failed for " + "; ".join(errors[:8]))

        for bam in (work / "out_dir").rglob("*.bam"):
            bam.unlink(missing_ok=True)

        tar_path = work / f"{sample_id}.locityper.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(work / "out_dir", arcname="out_dir")
        n_rows = write_summary(work, sample_id, work / "gts.filtered.csv")
        print(f"wrote {n_rows} loci to gts.filtered.csv", flush=True)
        print(f"wrote {tar_path} ({tar_path.stat().st_size} bytes)", flush=True)
        print("genotype done", flush=True)
    except Exception:
        status = 1
        traceback.print_exc()
    finally:
        _dump_monitor(monitor, work, sample_id=sample_id, stage="genotype")
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_fh.close()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
