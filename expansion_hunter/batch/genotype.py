#!/usr/bin/env python3
"""Native Batch genotype worker (EH 0.1.0: ExpansionHunter + python3, no GCS).

Inputs are already on /work (host runnable downloaded them). Writes
{sample}.EH.json and {sample}.EH.vcf next to the minicram. Host runnable
uploads those after this process exits.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import traceback
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


def sex_flag(raw: str) -> str:
    key = raw.lower().replace(" ", "")
    if key in {"male", "m", "1"}:
        return "male"
    if key in {"female", "f", "2"}:
        return "female"
    print(f"WARN: unrecognized sex {raw!r}; defaulting to female", flush=True)
    return "female"


def eh_help() -> str:
    proc = subprocess.run(
        ["ExpansionHunter", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    return (proc.stdout or "") + (proc.stderr or "")


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
        prefix = f"{sample_id}.EH"
        reads = work / "reads.cram"
        reads_idx = work / "reads.cram.crai"
        ref = work / "reference.fa"
        catalog = work / "catalog.json"
        for path in (reads, reads_idx, ref, work / "reference.fa.fai", catalog):
            if not path.is_file() or path.stat().st_size <= 0:
                raise SystemExit(f"missing input: {path}")

        json.load(catalog.open())
        sex = sex_flag(os.environ.get("SEX", "female"))
        help_text = eh_help()
        cat_flag = "--variant-catalog" if "--variant-catalog" in help_text else "--catalog"
        mode = "optimized-streaming" if "optimized-streaming" in help_text else "seeking"
        extra = [
            flag
            for flag in (
                "--dont-output-consensus-sequences",
                "--dont-output-quality-metrics",
                "--disable-all-plots",
            )
            if flag in help_text
        ]
        eh_bin = shutil.which("ExpansionHunter")
        if not eh_bin:
            raise SystemExit("ExpansionHunter is not on PATH")
        print(f"genotype start sample={sample_id} sex={sex} bin={eh_bin}", flush=True)
        subprocess.run(["ExpansionHunter", "--version"], check=False)
        print(f"catalog flag={cat_flag} analysis-mode={mode}", flush=True)

        cmd = [
            "ExpansionHunter",
            "--reads",
            str(reads),
            "--reads-index",
            str(reads_idx),
            "--reference",
            str(ref),
            cat_flag,
            str(catalog),
            "--output-prefix",
            prefix,
            "--sex",
            sex,
            "--analysis-mode",
            mode,
            *extra,
        ]
        print("+", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=work)

        out_json = work / f"{prefix}.json"
        out_vcf = work / f"{prefix}.vcf"
        if not out_json.is_file() or out_json.stat().st_size <= 0:
            raise SystemExit(f"ExpansionHunter wrote no {out_json}")
        if not out_vcf.is_file() or out_vcf.stat().st_size <= 0:
            raise SystemExit(f"ExpansionHunter wrote no {out_vcf}")
        json.load(out_json.open())
        print(f"wrote {out_json} ({out_json.stat().st_size} bytes)", flush=True)
        print(f"wrote {out_vcf} ({out_vcf.stat().st_size} bytes)", flush=True)
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
