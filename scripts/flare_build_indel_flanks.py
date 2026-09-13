#!/usr/bin/env python3
"""Build a BED of indel ± flank intervals from a (possibly remote) VCF.

Used when ``indel_flank_bp > 0`` so target and reference can exclude the same
padded indel neighborhoods. Emits 0-based half-open BED.

    python3 scripts/flare_build_indel_flanks.py \\
      --vcf gt.vcf.gz --region chr20 --flank 5 --out indel_flanks.bed
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional


def iter_indel_intervals(vcf: str, region: str, flank: int):
    # bcftools query has no -v type filter (-v is verbosity). Use -i TYPE=.
    cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\n",
        "-i",
        'TYPE="indel"',
    ]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(vcf)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
        pos = int(pos_s)
        # Span covering REF and all ALTs (VCF POS is 1-based start of REF).
        allele_lens = [len(ref)] + [len(a) for a in alt.split(",")]
        end_1based = pos + max(allele_lens) - 1
        start0 = max(0, pos - 1 - flank)
        end0 = end_1based + flank
        yield chrom, start0, end0
    stderr = proc.communicate()[1]
    rc = proc.returncode
    if rc != 0:
        if stderr:
            print(stderr, file=sys.stderr, end="")
        raise SystemExit(f"bcftools query failed ({rc})")


def write_bed(rows, dest: Path) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with dest.open("w") as fh:
        for chrom, start, end in rows:
            if end <= start:
                continue
            fh.write(f"{chrom}\t{start}\t{end}\n")
            n += 1
    return n


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--region", default="")
    p.add_argument("--flank", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)
    if args.flank < 0:
        raise SystemExit("--flank must be >= 0")
    if args.flank == 0:
        args.out.write_text("")
        print("flank=0; empty bed", args.out)
        return 0
    n = write_bed(iter_indel_intervals(args.vcf, args.region, args.flank), args.out)
    print(f"wrote {n} indel-flank intervals -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
