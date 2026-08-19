#!/usr/bin/env python3
"""
Fill AC/AN/AF and per-sample carrier presence from a multi-sample VCF.

No-call policy: missing genotypes do not contribute alleles and are not carriers.
Writes:
  - updated site table (freq fields)
  - carriers TSV: site_id \\t sample1,sample2,...

Streams the site table in VCF order (same order as extract_sites_from_vcf.py)
so a multi-million-site callset does not have to sit in RAM.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import (  # noqa: E402
    freq_class,
    iter_gt_records,
    iter_sites,
    list_vcf_samples,
    open_text,
    site_row,
    write_site_header,
)

# Common diploid / missing GTs — avoid split/parse on every sample of every site.
_GT_FAST: dict[str, tuple[int, int, bool]] = {
    "0/0": (0, 2, False),
    "0|0": (0, 2, False),
    "0/1": (1, 2, True),
    "1/0": (1, 2, True),
    "0|1": (1, 2, True),
    "1|0": (1, 2, True),
    "1/1": (2, 2, True),
    "1|1": (2, 2, True),
    "./.": (0, 0, False),
    ".|.": (0, 0, False),
    ".": (0, 0, False),
    "0": (0, 1, False),
    "1": (1, 1, True),
}


def gt_stats(gt: str) -> tuple[int, int, bool]:
    """Return (alt_count, allele_number, is_carrier)."""
    hit = _GT_FAST.get(gt)
    if hit is not None:
        return hit
    if not gt or gt.startswith("./.") or gt.startswith(".|."):
        return 0, 0, False
    alleles = gt.replace("|", "/").split("/")
    if all(a == "." for a in alleles):
        return 0, 0, False
    called = [a for a in alleles if a != "."]
    if not called:
        return 0, 0, False
    alt = sum(1 for a in called if a != "0")
    return alt, len(called), alt > 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="Multi-sample VCF/VCF.gz/BCF")
    p.add_argument("--sites", required=True, help="Site table from extract_sites_from_vcf.py")
    p.add_argument("--out-sites", required=True)
    p.add_argument("--out-carriers", required=True)
    args = p.parse_args()

    sample_names = list_vcf_samples(args.vcf)
    n_samples = len(sample_names)
    print(f"[fill_af] {n_samples:,} samples; streaming in VCF order", file=sys.stderr, flush=True)

    Path(args.out_sites).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_carriers).parent.mkdir(parents=True, exist_ok=True)

    sites_iter = iter_sites(args.sites)
    n = 0
    with open_text(args.out_sites, "wt") as out_s, open_text(args.out_carriers, "wt") as out_c:
        write_site_header(out_s)
        out_c.write("id\tcarriers\tn_carriers\n")
        for chrom, pos, vid, gts in iter_gt_records(args.vcf):
            try:
                site = next(sites_iter)
            except StopIteration as exc:
                raise RuntimeError(
                    f"VCF has more records than the sites TSV at {chrom}:{pos} {vid}"
                ) from exc
            if str(site["chrom"]) != chrom or int(site["pos"]) != pos:
                raise RuntimeError(
                    f"site/VCF order mismatch at record {n + 1}: "
                    f"site {site['chrom']}:{site['pos']} {site['id']} vs "
                    f"VCF {chrom}:{pos} {vid}"
                )
            n += 1
            if n % 50000 == 0:
                print(f"[fill_af] {n:,} records", file=sys.stderr, flush=True)

            site_ac = 0
            site_an = 0
            carriers: list[str] = []
            limit = min(len(gts), n_samples)
            for si in range(limit):
                ac_i, an_i, is_car = gt_stats(gts[si])
                site_ac += ac_i
                site_an += an_i
                if is_car:
                    carriers.append(sample_names[si])

            site["ac"] = site_ac
            site["an"] = site_an
            site["af"] = (site_ac / site_an) if site_an else 0.0
            site["n_carriers"] = len(carriers)
            site["freq_class"] = freq_class(site_ac, site_an, site["n_carriers"], n_samples)
            out_s.write(site_row(site))
            out_c.write(f"{site['id']}\t{','.join(carriers)}\t{site['n_carriers']}\n")

        extra = next(sites_iter, None)
        if extra is not None:
            raise RuntimeError(
                f"sites TSV has extra rows after the VCF ended, starting at {extra['id']}"
            )

    print(f"[fill_af] filled {n:,} records", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
