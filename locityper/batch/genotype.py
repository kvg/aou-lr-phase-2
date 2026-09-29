#!/usr/bin/env python3
"""Host-side Locityper summary after the container worker writes out_dir.

The Locityper image has locityper + bash, not python3 (WDL Summarize uses
util_docker). Batch runs this on the VM host with LT_WORKDIR=/mnt/disks/lt.
"""

from __future__ import annotations

import gzip
import json
import math
import os
from pathlib import Path


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
    work = Path(os.environ.get("LT_WORKDIR", "/work"))
    sample_id = os.environ.get("SAMPLE_ID", "").strip()
    if not sample_id:
        raise SystemExit("SAMPLE_ID is empty")
    tar_path = work / f"{sample_id}.locityper.tar.gz"
    if not tar_path.is_file() or tar_path.stat().st_size <= 0:
        raise SystemExit(f"missing {tar_path}")
    n_rows = write_summary(work, sample_id, work / "gts.filtered.csv")
    print(f"wrote {n_rows} loci to gts.filtered.csv", flush=True)
    print(f"tar {tar_path} ({tar_path.stat().st_size} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
