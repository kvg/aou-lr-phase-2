#!/usr/bin/env python3
"""Build the Terra entity table ``TrgtHaplotypePhase`` runs on.

The workflow needs four per-sample paths: the haplotagged BAM and its index,
and the TRGT VCF and its index. They live in two different places — the BAMs in
a Terra table (``aou2_v1_phased_bams``), the TRGT calls in
``tractor_mix/resources/trgt_table.txt`` — and neither is a superset of the
other. This joins them and emits one row per sample that has **both**.

    make_trgt_phase_table.py --bams-table aou2_v1_phased_bams.tsv \\
        --trgt-table tractor_mix/resources/trgt_table.txt \\
        --out trgt_phase.tsv

Everything else the workflow takes is the same for every sample — the phased
per-chromosome callsets are multi-sample files, and the reference, catalog,
thresholds and image are run-level — so they belong in the method config, not
in this table. A column that is identical in 12,000 rows is a column that can
disagree with itself.

Rows are emitted only for the intersection, and the three populations are
reported separately: a sample with a BAM but no TRGT call cannot be phased, and
a TRGT call with no BAM has no reads to phase it from. Neither is an error, but
both are a count you want to see rather than discover as a failed shard.

``--on-duplicate`` is required when an emitted sample has more than one TRGT
path. ``trgt_table.txt`` holds 12,421 rows but 12,393 distinct person ids: 28
samples were genotyped under two different Terra submissions and both copies are
listed. One column holds one URI, so the choice is explicit rather than decided
by file order.

Nothing here checks that the objects exist. It is a join over path strings.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Optional, TextIO

OUT_COLUMNS = ("sample", "bam", "bam_index", "trgt_vcf", "trgt_vcf_index")


def open_text(path: str, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_trgt_table(path: str) -> tuple[dict[str, list[str]], list[str]]:
    """(person_id -> VCF URIs in file order, ids with more than one distinct URI).

    Reading never fails on a duplicate. Whether a duplicate matters depends on
    whether that id is actually emitted, which only the caller knows, so the
    policy is applied in ``resolve`` and enforced in ``run``.
    """
    out: dict[str, list[str]] = {}
    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise SystemExit(f"{path}:{n}: expected 'id<TAB>uri', got {line!r}")
            sample, uri = fields[0].strip(), fields[1].strip()
            paths = out.setdefault(sample, [])
            if uri not in paths:
                paths.append(uri)
    dups = [s for s, paths in out.items() if len(paths) > 1]
    return out, dups


def resolve(paths: dict[str, list[str]], on_duplicate: str) -> dict[str, str]:
    """One URI per id. ``fail`` resolves like ``first``; ``run`` refuses before it is used."""
    idx = -1 if on_duplicate == "last" else 0
    return {sample: uris[idx] for sample, uris in paths.items()}


def read_bams_table(path: str, bam_column: str, bai_column: str) -> dict[str, dict[str, str]]:
    """entity id -> {bam, bam_index}, from a Terra table export.

    The id is the first column, whose header Terra writes as ``entity:<type>_id``.
    Rows whose BAM or index cell is empty are dropped and counted by the caller:
    an empty cell in Terra is a real state, not a missing file.
    """
    rows: dict[str, dict[str, str]] = {}
    with open_text(path) as fh:
        header: Optional[list[str]] = None
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = [f.strip() for f in fields]
                if not header[0].startswith("entity:"):
                    raise SystemExit(
                        f"{path}: first column header is {header[0]!r}, expected 'entity:<type>_id'. "
                        "Export the table from Terra rather than hand-building it.")
                for col in (bam_column, bai_column):
                    if col not in header:
                        raise SystemExit(f"{path}: no column {col!r}. Columns: {', '.join(header)}")
                continue
            row = dict(zip(header, fields))
            sample = fields[0].strip()
            if not sample:
                continue
            rows[sample] = {
                "bam": row.get(bam_column, "").strip(),
                "bam_index": row.get(bai_column, "").strip(),
            }
    if header is None:
        raise SystemExit(f"{path}: empty")
    return rows


def run(args: argparse.Namespace) -> int:
    bams = read_bams_table(args.bams_table, args.bam_column, args.bai_column)
    paths, dups = read_trgt_table(args.trgt_table)
    trgt = resolve(paths, args.on_duplicate)

    incomplete = sorted(s for s, r in bams.items() if not r["bam"] or not r["bam_index"])
    usable_bams = {s: r for s, r in bams.items() if r["bam"] and r["bam_index"]}

    both = [s for s in usable_bams if s in trgt]
    bam_only = sorted(s for s in usable_bams if s not in trgt)
    trgt_only = sorted(s for s in trgt if s not in usable_bams)

    conflicting = [s for s in both if s in dups]
    if conflicting and args.on_duplicate == "fail":
        s = conflicting[0]
        listing = "\n".join(f"  {u}" for u in paths[s])
        raise SystemExit(
            f"{len(conflicting)} emitted sample(s) carry more than one TRGT path, e.g. {s}:\n"
            f"{listing}\nall: {', '.join(conflicting)}\n"
            "Choose with --on-duplicate first|last.")

    # Emit in the BAM table's order so the new table lines up with the old one.
    order = [s for s in bams if s in set(both)]
    with open(args.out, "wt") as fh:
        fh.write(f"entity:{args.entity_type}_id\t" + "\t".join(OUT_COLUMNS) + "\n")
        for sample in order:
            fh.write("\t".join((
                sample, sample, usable_bams[sample]["bam"], usable_bams[sample]["bam_index"],
                trgt[sample], trgt[sample] + args.index_suffix)) + "\n")

    print(f"BAM table rows:                 {len(bams):,}")
    print(f"  with both BAM and index:      {len(usable_bams):,}")
    print(f"  missing a BAM or index cell:  {len(incomplete):,}"
          + (f"  (e.g. {', '.join(incomplete[:5])})" if incomplete else ""))
    print(f"distinct TRGT samples:          {len(trgt):,}")
    print(f"  ids with >1 TRGT path:        {len(dups):,}")
    print()
    print(f"ROWS WRITTEN (both sources):    {len(order):,}")
    print(f"  BAM but no TRGT call:         {len(bam_only):,}"
          + (f"  (e.g. {', '.join(bam_only[:5])})" if bam_only else ""))
    print(f"  TRGT call but no usable BAM:  {len(trgt_only):,}"
          + (f"  (e.g. {', '.join(trgt_only[:5])})" if trgt_only else ""))
    if conflicting:
        print(f"\nresolved by --on-duplicate {args.on_duplicate}: "
              f"{', '.join(conflicting[:10])}" + (" ..." if len(conflicting) > 10 else ""),
              file=sys.stderr)
    print(f"\nwrote {args.out}")
    if bam_only or trgt_only:
        print("Samples outside the intersection are absent from this table by construction, so "
              "they cannot be submitted by accident.", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bams-table", required=True,
                   help="Terra export of the haplotagged-BAM table (aou2_v1_phased_bams)")
    p.add_argument("--trgt-table", required=True,
                   help="person_id<TAB>TRGT VCF URI, no header")
    p.add_argument("--entity-type", default="trgt_phase",
                   help="Entity type for the new table, used in the 'entity:<type>_id' header")
    p.add_argument("--bam-column", default="gcs_bam_path")
    p.add_argument("--bai-column", default="gcs_bai_path")
    p.add_argument("--index-suffix", default=".tbi")
    p.add_argument("--on-duplicate", choices=("fail", "first", "last"), default="fail")
    p.add_argument("--out", required=True)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
