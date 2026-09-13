#!/usr/bin/env python3
"""Extract site-level rows from an SV VCF partition (no genotypes required)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running as scripts/foo.py
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import (  # noqa: E402
    abs_svlen,
    cadd_sv_bin,
    infer_svlen,
    infer_svtype,
    iter_site_records,
    open_text,
    parse_info,
    site_row,
    size_flags,
    write_site_header,
)


def compact_alt(alt: str, svtype: str) -> str:
    if alt and alt != "." and len(alt) <= 64:
        return alt
    if svtype:
        return f"<{svtype}>"
    return alt if alt else "."


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--source-vcf", required=True, choices=["main", "bnd", "large"])
    p.add_argument("--phase", required=True, choices=["phase1", "phase2"])
    p.add_argument("--out", required=True)
    args = p.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    with open_text(out_path, "wt") as out:
        write_site_header(out)
        for chrom, pos, vid, ref, alt, filt, info in iter_site_records(args.vcf):
            n += 1
            if n % 50000 == 0:
                print(f"[extract_sites] {n:,} records", file=sys.stderr, flush=True)
            info_d = parse_info(info)
            svtype = infer_svtype(alt, info_d, args.source_vcf)
            svlen = infer_svlen(ref, alt, info_d)
            if svlen is None and "END" in info_d:
                try:
                    end = int(info_d["END"])
                    svlen = end - int(pos)
                    if svtype == "DEL":
                        svlen = -abs(svlen)
                except ValueError:
                    end = None
            else:
                end = None
            if end is None and "END" in info_d:
                try:
                    end = int(info_d["END"])
                except ValueError:
                    end = int(pos)
            if end is None:
                end = int(pos) + abs(svlen or 0)

            ge20, ge50 = size_flags(svlen, svtype)
            row = {
                "chrom": chrom,
                "pos": int(pos),
                "end": end,
                "id": vid if vid != "." else f"{chrom}:{pos}:{svtype}",
                "ref": ref,
                "alt": compact_alt(alt, svtype),
                "svtype": svtype,
                "svlen": svlen if svlen is not None else "",
                "filter": filt,
                "source_vcf": args.source_vcf,
                "phase": args.phase,
                "ac": 0,
                "an": 0,
                "af": 0.0,
                "n_carriers": 0,
                "freq_class": "singleton",
                "region_class": "non_repetitive",
                "hit_rmsk": False,
                "hit_simpleRepeat": False,
                "hit_genomicSuperDups": False,
                "hit_cmrg": False,
                "cadd_sv_phred": "",
                "cadd_sv_bin": cadd_sv_bin(None),
                "size_bin_ge20": ge20,
                "size_bin_ge50": ge50,
                "suppressed_main_duplicate": False,
            }
            _ = abs_svlen(svlen)
            out.write(site_row(row))
    print(f"[extract_sites] wrote {n:,} sites to {out_path}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
