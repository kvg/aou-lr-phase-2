#!/usr/bin/env python3
"""Rebuild TRGT's repeat-definition BED from a TRGT output VCF.

``trgt --repeats`` needs four tab-separated fields::

    chrom  start  end  ID=<trid>;MOTIFS=<motifs>;STRUC=<struc>

The catalog in ``refs/repeat_loci/`` is five fields (``chrom start end TRID
motifs``) because ``trgt_plvi.py catalog`` wrote it for
``aggregate_repeat_loci.py``, dropping ``STRUC`` on the way. ``STRUC`` cannot be
reconstructed from a motif list at a compound locus, so that file cannot be
converted — but a TRGT *output* VCF carries all three fields in INFO, because
TRGT echoes the definition it was given::

    TRID=1-66431-66446-TATAT;END=66446;MOTIFS=TATAT;STRUC=<TR22>
    TRID=1-66160-66175-TA,1-66192-66201-TA,...;END=66415;MOTIFS=ATATA,TATAT,TA;STRUC=<VC2>

So one genotyped sample round-trips the definitions it was genotyped against,
which is also the set that makes a re-genotyping comparison locus-for-locus.

    trgt_defs_from_vcf.py --vcf <sample>_trgt.TR_Explorer_1.01.vcf.gz \\
        --out trexplorer_v1.0.1.trgt_definitions.bed

Coordinates: VCF ``POS`` is the 1-based locus start and BED is half-open, so the
start is ``POS - 1`` and the end is INFO ``END``. This matches the derived
catalog, whose rows read ``chr20 10000929 10000950 20-10000930-10000950-A``.

**One sample's VCF holds the loci TRGT emitted for that sample.** Any locus it
dropped is absent here. Check the row count against the catalog ceiling
(``gsutil cat <catalog> | gunzip | wc -l``) before treating the output as the
whole catalog; ``--expect-loci`` turns that check into a hard one.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path
from typing import Iterator, Optional

QUERY = "%CHROM\t%POS\t%INFO/END\t%INFO/TRID\t%INFO/MOTIFS\t%INFO/STRUC\n"
MISSING = (".", "")


def query_lines(vcf: str) -> Iterator[str]:
    """bcftools query stream. stderr is surfaced: an unreadable URI must not look empty."""
    proc = subprocess.Popen(["bcftools", "query", "-f", QUERY, vcf],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            bufsize=1 << 20)
    assert proc.stdout is not None and proc.stderr is not None
    try:
        yield from proc.stdout
    finally:
        proc.stdout.close()
        err = proc.stderr.read()
        proc.stderr.close()
        if proc.wait() != 0:
            raise SystemExit(f"bcftools query failed on {vcf}:\n{err.strip()}")


def definition(chrom: str, pos: str, end: str, trid: str, motifs: str, struc: str) -> str:
    return f"{chrom}\t{int(pos) - 1}\t{int(end)}\tID={trid};MOTIFS={motifs};STRUC={struc}"


def run(args: argparse.Namespace) -> dict:
    counts = {"records": 0, "written": 0, "skipped_missing_field": 0}
    skipped: list[str] = []
    opener = gzip.open if args.out.endswith(".gz") else open
    with opener(args.out, "wt") as fh:
        for line in query_lines(args.vcf):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            counts["records"] += 1
            chrom, pos, end, trid, motifs, struc = fields[:6]
            if any(v in MISSING for v in (end, trid, motifs, struc)):
                counts["skipped_missing_field"] += 1
                if len(skipped) < 5:
                    skipped.append(f"{chrom}:{pos} TRID={trid} MOTIFS={motifs} STRUC={struc}")
                continue
            fh.write(definition(chrom, pos, end, trid, motifs, struc) + "\n")
            counts["written"] += 1

    if counts["written"] == 0:
        raise SystemExit(
            f"error: no definition written from {args.vcf}. If every record lacked STRUC, this "
            "is not a TRGT output VCF, or TRGT was run without a STRUC-bearing catalog.")
    if counts["skipped_missing_field"]:
        print(f"warning: skipped {counts['skipped_missing_field']:,} record(s) missing END, TRID, "
              f"MOTIFS or STRUC; e.g. {'; '.join(skipped)}", file=sys.stderr)
    if args.expect_loci is not None and counts["written"] != args.expect_loci:
        raise SystemExit(
            f"error: wrote {counts['written']:,} definitions but --expect-loci is "
            f"{args.expect_loci:,}. One sample's VCF holds only the loci TRGT emitted for it.")
    return counts


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="A TRGT output VCF (local path or gs:// URI)")
    p.add_argument("--out", required=True, help="Output BED; .gz is compressed")
    p.add_argument("--expect-loci", type=int, default=None,
                   help="Fail unless exactly this many definitions are written, e.g. the "
                        "catalog's line count")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    counts = run(args)
    for k, v in counts.items():
        print(f"{k:>24}  {v:,}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
