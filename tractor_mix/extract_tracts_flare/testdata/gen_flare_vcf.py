#!/usr/bin/env python3
"""Generate a deterministic FLARE-like VCF for parity and benchmark tests."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def open_out(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return path.open("w")


def allele(site: int, sample: int, hap: int) -> str:
    return str((site + sample + hap) % 2)


def ancestry(site: int, sample: int, hap: int, num_ancs: int) -> str:
    return str((site * 3 + sample * 5 + hap * 7) % num_ancs)


def generate(path: Path, n_sites: int, n_samples: int, num_ancs: int, with_anp: bool) -> None:
    samples = [f"S{i}" for i in range(n_samples)]
    fmt = "GT:AN1:AN2:ANP1:ANP2" if with_anp else "GT:AN1:AN2"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_out(path) as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry hap 1">\n')
        fh.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry hap 2">\n')
        for i in range(num_ancs):
            fh.write(f"##ANCESTRY=<ID=anc{i},Value={i}>\n")
        fh.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        for site in range(n_sites):
            pos = 1000 + site
            fields = [
                "chr22",
                str(pos),
                f"rs{site}",
                "A",
                "G",
                ".",
                "PASS",
                ".",
                fmt,
            ]
            for s in range(n_samples):
                a = allele(site, s, 0)
                b = allele(site, s, 1)
                c0 = ancestry(site, s, 0, num_ancs)
                c1 = ancestry(site, s, 1, num_ancs)
                rec = f"{a}|{b}:{c0}:{c1}"
                if with_anp:
                    rec += ":0.8,0.2:0.7,0.3"
                fields.append(rec)
            fh.write("\t".join(fields) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--sites", type=int, default=20)
    p.add_argument("--samples", type=int, default=8)
    p.add_argument("--num-ancs", type=int, default=5)
    p.add_argument("--anp", action="store_true")
    args = p.parse_args()
    generate(args.out, args.sites, args.samples, args.num_ancs, args.anp)
    print(f"Wrote {args.out} ({args.sites} sites x {args.samples} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
