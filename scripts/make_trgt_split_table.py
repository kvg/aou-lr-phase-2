#!/usr/bin/env python3
"""Terra table and sample-batches file for SubsetAndSplitVcf, per chromosome.

``TrgtHaplotypePhase`` currently extracts one sample from the multi-sample
ligated callset inside every task. Genome-wide over the cohort that is
12,393 x 22 = 272,646 passes over the phased callset. Splitting it once into
per-sample VCFs makes it 22, and every task then reads a few MB.

This writes the inputs for that split:

    --out-table    one row per autosome for SubsetAndSplitVcf
    --out-batches  the sample-batches file, so +split writes a few hundred
                   files per task instead of 12,393

Shape, and why it differs from split_short_regions
--------------------------------------------------
The earlier table sharded each chromosome 32 ways on site count, against a
joint-calling gVCF. Do not copy that here. A shard writes a per-sample file for
*its region*, so sub-chromosomal shards would give 32 x 22 x 12,393 = 8.7M
files and force a per-sample concatenation before anything could read them. One
row per chromosome gives each sample exactly one file per chromosome, which is
what the scatter wants.

The scale problem moves to the sample axis instead, and ``bcftools +split``
handles it through ``--samples-file``: the two-stage path in SubsetAndSplitVcf
splits into batches first, then each batch into per-sample files. That divides
the file count per task without dividing the region.

Sample set
----------
``--samples`` must be the same cohort the run uses, so pass the ``trgt_phase``
table (or its id column). Splitting samples that will never be genotyped wastes
storage; missing one means a task fails on a file that was never written.

Outputs are named ``<sample>.<output_tag>.vcf.gz`` by the workflow, with
``output_tag = <prefix>.<chrom>``, so a sample's chromosome file is at
``<out-base>/<chrom>/<sample>.<prefix>.<chrom>.vcf.gz`` -- deterministic, so the
trgt_phase array columns can be generated without listing the bucket.

bgzipped VCF rather than BCF because whatshap reads these files directly and
VCF is the format it documents; the workflow indexes every one, since a region
read and ``whatshap --regions`` both need random access.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Optional, TextIO

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
TABLE_COLUMNS = ("region", "joint_vcf", "joint_vcf_idx", "gcs_output_dir",
                 "output_tag", "sample_batches_tsv")


def open_text(path: str, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_lai_table(path: str, chrom_col: str, vcf_col: str, idx_col: str
                   ) -> dict[str, tuple[str, str]]:
    """chrom -> (vcf, index), from a flare_lai_prod export."""
    out: dict[str, tuple[str, str]] = {}
    with open_text(path) as fh:
        header: Optional[list[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = [f.strip() for f in fields]
                missing = [c for c in (chrom_col, vcf_col, idx_col) if c not in header]
                if missing:
                    raise SystemExit(
                        f"{path}: no column(s) {', '.join(missing)}. Columns: {', '.join(header)}")
                continue
            row = dict(zip(header, fields))
            chrom = row.get(chrom_col, "").strip()
            vcf, idx = row.get(vcf_col, "").strip(), row.get(idx_col, "").strip()
            if not chrom or not vcf or not idx:
                continue
            out[chrom] = (vcf, idx)
    if header is None:
        raise SystemExit(f"{path}: empty")
    return out


def read_samples(path: str) -> list[str]:
    """Entity ids in file order, from an id list or a Terra table export."""
    ids, seen = [], set()
    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            first = line.split("\t")[0].strip()
            if n == 1 and first.startswith("entity:"):
                continue
            if first and first not in seen:
                seen.add(first)
                ids.append(first)
    return ids


def batch_rows(samples: list[str], batch_size: int) -> list[tuple[str, str]]:
    """(batch name, comma-joined samples) in blocks of ``batch_size``.

    Names are ``batch-<i>`` with NO zero padding, because SubsetAndSplitVcf's
    scatter reconstructs each path as ``"batch-" + i`` from a plain Int. Padding
    them to ``batch-00`` broke exactly the first ten: 10..24 happen to agree,
    0..9 do not, and the workflow failed with FileNotFoundException on precisely
    those. The lexical-sort tidiness padding buys is worth nothing here -- the
    workflow indexes by integer and never sorts the names.
    """
    out = []
    for i in range(0, len(samples), batch_size):
        block = samples[i:i + batch_size]
        out.append((f"batch-{i // batch_size}", ",".join(block)))
    return out


def run(args: argparse.Namespace) -> int:
    lai = read_lai_table(args.lai_table, args.chrom_column, args.vcf_column, args.index_column)
    samples = read_samples(args.samples)
    if not samples:
        raise SystemExit(f"{args.samples}: no sample ids")

    chroms = [c for c in AUTOSOMES if c in lai]
    missing = [c for c in AUTOSOMES if c not in lai]
    extra = sorted(set(lai) - set(AUTOSOMES))

    # The sample-batches file. Format is the one SubsetAndSplitVcf documents:
    # comma-joined samples, a '-' placeholder, then the batch name.
    batches = batch_rows(samples, args.batch_size)
    with open(args.out_batches, "wt") as fh:
        for name, joined in batches:
            fh.write(f"{joined}\t-\t{name}\n")

    base = args.out_base.rstrip("/")
    with open(args.out_table, "wt") as fh:
        fh.write(f"entity:{args.entity_type}_id\t" + "\t".join(TABLE_COLUMNS) + "\n")
        for chrom in chroms:
            vcf, idx = lai[chrom]
            fh.write("\t".join((
                chrom, chrom, vcf, idx, f"{base}/{chrom}",
                f"{args.output_prefix}.{chrom}", args.batches_uri or args.out_batches)) + "\n")

    print(f"samples:            {len(samples):,}")
    print(f"batches:            {len(batches):,} of up to {args.batch_size:,}")
    print(f"files per task:     {min(args.batch_size, len(samples)):,} "
          f"(vs {len(samples):,} unbatched)")
    print(f"table rows:         {len(chroms)}")
    if missing:
        print(f"autosomes absent from the LAI table: {', '.join(missing)}", file=sys.stderr)
    if extra:
        print(f"non-autosomes ignored: {', '.join(extra)}", file=sys.stderr)
    print(f"\nper-sample file path will be:\n  {base}/<chrom>/<sample>."
          f"{args.output_prefix}.<chrom>.vcf.gz")
    print(f"\nwrote {args.out_table}\nwrote {args.out_batches}")
    if not args.batches_uri:
        print("\nUpload the batches file and re-run with --batches-uri, or edit the "
              "sample_batches_tsv column: Terra needs a gs:// path, not a local one.",
              file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lai-table", required=True,
                   help="flare_lai_prod export carrying chrom, gt_vcf, gt_vcf_index")
    p.add_argument("--samples", required=True,
                   help="trgt_phase table export, or a bare id list")
    p.add_argument("--out-base", required=True,
                   help="gs:// prefix for the per-sample VCFs; <out-base>/<chrom>/ per row")
    p.add_argument("--batch-size", type=int, default=500,
                   help="Samples per +split batch (default 500)")
    p.add_argument("--output-prefix", default="phased",
                   help="output_tag is <prefix>.<chrom> (default phased)")
    p.add_argument("--batches-uri", default=None,
                   help="gs:// path the batches file will live at; used in the table column")
    p.add_argument("--entity-type", default="trgt_split_regions")
    p.add_argument("--chrom-column", default="chrom")
    p.add_argument("--vcf-column", default="gt_vcf")
    p.add_argument("--index-column", default="gt_vcf_index")
    p.add_argument("--out-table", required=True)
    p.add_argument("--out-batches", required=True)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
