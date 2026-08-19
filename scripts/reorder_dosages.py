#!/usr/bin/env python3
"""Subset and reorder Tractor dosage/hapcount columns to match analysis_samples.txt."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_ids(path: Path) -> list[str]:
    ids = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s:
                ids.append(s.split()[0])
    return ids


def reorder_file(infile: Path, outfile: Path, sample_order: list[str]) -> None:
    with open_text(infile, "rt") as fin:
        header = fin.readline().rstrip("\n").split("\t")
        if len(header) < 6:
            raise SystemExit(f"{infile}: expected meta + sample columns, got {header[:10]}")
        meta = header[:5]
        samples = header[5:]
        index = {s: i for i, s in enumerate(samples)}
        missing = [s for s in sample_order if s not in index]
        if missing:
            raise SystemExit(
                f"{infile}: {len(missing)} analysis samples missing from dosage "
                f"(e.g. {missing[:5]})"
            )
        cols = [index[s] for s in sample_order]
        with open_text(outfile, "wt") as fout:
            fout.write("\t".join(meta + sample_order) + "\n")
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                row_meta = parts[:5]
                vals = parts[5:]
                fout.write("\t".join(row_meta + [vals[i] for i in cols]) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True, type=Path)
    p.add_argument(
        "--in-files",
        nargs="+",
        required=True,
        type=Path,
        help="Dosage and/or hapcount files to reorder",
    )
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()

    samples = read_ids(args.samples)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for inf in args.in_files:
        # Keep .gz suffix if present
        out = args.out_dir / inf.name.replace(".dosage.txt", ".ordered.dosage.txt").replace(
            ".hapcount.txt", ".ordered.hapcount.txt"
        )
        if out == args.out_dir / inf.name:
            out = args.out_dir / f"ordered.{inf.name}"
        print(f"Reordering {inf} -> {out} ({len(samples)} samples)")
        reorder_file(inf, out, samples)
    # Write checksum of sample order for downstream asserts
    (args.out_dir / "dosage_sample_order.txt").write_text("\n".join(samples) + "\n")


if __name__ == "__main__":
    main()
