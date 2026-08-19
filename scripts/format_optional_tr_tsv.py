#!/usr/bin/env python3
"""
Optional: convert AnnotateIndelTRs-style TSV into annotate-ready columns.

Input TSV columns (no header): CHROM POS REF ALT ID DETECTED MOTIF END START_0BASED MOTIF_SIZE
Output: CHROM POS REF ALT ID TR_DETECTED TR_MOTIF TR_END TR_START_0BASED TR_MOTIF_SIZE
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-tsv", required=True)
    p.add_argument("--out-tsv", required=True)
    args = p.parse_args()

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.in_tsv) as inf, open(args.out_tsv, "w") as out:
        out.write(
            "CHROM\tPOS\tREF\tALT\tID\tTR_DETECTED\tTR_MOTIF\tTR_END\tTR_START_0BASED\tTR_MOTIF_SIZE\n"
        )
        for line in inf:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            # pad optional motif fields
            while len(parts) < 10:
                parts.append(".")
            chrom, pos, ref, alt, vid = parts[:5]
            detected, motif, end, start0, msize = parts[5:10]
            out.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{vid}\t{detected}\t{motif}\t{end}\t{start0}\t{msize}\n"
            )


if __name__ == "__main__":
    main()
