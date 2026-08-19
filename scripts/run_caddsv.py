#!/usr/bin/env python3
"""Run CADD-SV v2 with a staged annotation bundle and expose its score TSV."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

# CADD-SV prep_files calls pysam.fetch() on the full DUP/INV span before it
# truncates to ~1 kb. Ultralong intervals OOM (exit 137). DELs only need flanks.
MAX_DUP_INV_BP = 1_000_000


def extract_annotations(bundle: Path, destination: Path) -> None:
    """Extract the annotation directory, refusing archive traversal."""
    with tarfile.open(bundle, "r:*") as archive:
        root = destination.resolve()
        for member in archive.getmembers():
            target = (destination / member.name).resolve()
            if target != root and root not in target.parents:
                raise ValueError(f"Unsafe path in CADD-SV annotation bundle: {member.name}")
        archive.extractall(destination, filter="data")


def filter_caddsv_bed(src: Path, dst: Path, max_dup_inv_bp: int = MAX_DUP_INV_BP) -> int:
    """Write intervals CADD-SV can sequence-fetch without loading megabase DUP/INV."""
    n_in = n_out = n_skip = 0
    with src.open("rt", encoding="utf-8") as inf, dst.open("wt", encoding="utf-8") as out:
        for line in inf:
            if not line.strip() or line.startswith("#"):
                continue
            n_in += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                span = int(parts[2]) - int(parts[1])
            except ValueError:
                continue
            if parts[3] in {"DUP", "INV"} and span > max_dup_inv_bp:
                n_skip += 1
                continue
            out.write(line if line.endswith("\n") else line + "\n")
            n_out += 1
    print(
        f"[run_caddsv] kept {n_out:,}/{n_in:,} intervals; "
        f"skipped {n_skip:,} DUP/INV > {max_dup_inv_bp:,} bp",
        file=sys.stderr,
        flush=True,
    )
    return n_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--annotations-tar", required=True)
    parser.add_argument("--out-scores", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--conda-prefix", default="/caddsv-conda")
    parser.add_argument("--max-dup-inv-bp", type=int, default=MAX_DUP_INV_BP)
    args = parser.parse_args()

    capped = Path("caddsv.capped.bed")
    n_kept = filter_caddsv_bed(Path(args.input_bed), capped, args.max_dup_inv_bp)
    if n_kept == 0:
        Path(args.out_scores).write_text("chrom\tstart\tend\ttype\tCADD-SV_PHRED\n")
        return

    annotations_root = Path("caddsv_annotations")
    extract_annotations(Path(args.annotations_tar), annotations_root)

    # Accept either an archive that contains annotations/ or its contents directly.
    annotations_dir = annotations_root / "annotations"
    if not annotations_dir.exists():
        annotations_dir = annotations_root

    output_dir = Path("caddsv_results")
    command = [
        "caddsv",
        "run",
        str(capped.resolve()),
        "--annotations-dir",
        str(annotations_dir),
        "--output-dir",
        str(output_dir),
        "--threads",
        str(args.threads),
        "--conda-prefix",
        args.conda_prefix,
    ]
    subprocess.run(command, check=True)

    scored = sorted((output_dir / "scored").glob("*_score.tsv"))
    if len(scored) != 1:
        raise RuntimeError(f"Expected one CADD-SV score TSV, found {scored}")
    shutil.copyfile(scored[0], args.out_scores)


if __name__ == "__main__":
    main()
