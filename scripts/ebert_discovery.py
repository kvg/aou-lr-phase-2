#!/usr/bin/env python3
"""
Build Ebert-style cumulative discovery tables.

Outputs a long TSV with cumulative counts after each ancestry-ordered sample,
optionally stratified by region_class and/or cadd_sv_bin.

Streams ``--sites`` and ``--carriers`` in lockstep (same order as
fill_af_and_carriers.py). Does not keep a per-site ID map in memory.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text, order_samples_by_ancestry  # noqa: E402


def load_ancestry(path: str) -> dict[str, str]:
    """Two-column TSV: sample_id \\t ancestry_label (header optional)."""
    out: dict[str, str] = {}
    with open_text(path, "rt") as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0].lower() in ("sample", "sample_id", "iid", "id"):
                continue
            out[parts[0]] = parts[1] if len(parts) > 1 else "oth"
    return out


def first_carrier_rank(carriers_field: str, sample_rank: dict[str, int]) -> int | None:
    if not carriers_field:
        return None
    best: int | None = None
    for sample in carriers_field.split(","):
        rank = sample_rank.get(sample)
        if rank is None:
            continue
        if best is None or rank < best:
            best = rank
    return best


def stratum_key(strata: str, region_class: str, cadd_sv_bin: str) -> str:
    if strata == "none":
        return "all"
    if strata == "region":
        return region_class
    if strata == "cadd":
        return cadd_sv_bin
    return f"{region_class}|{cadd_sv_bin}"


def write_discovery(
    out_path: Path,
    order: list[str],
    ancestry: dict[str, str],
    buckets: dict[tuple[str, str], list[int]],
) -> None:
    fieldnames = [
        "sample_index",
        "sample_id",
        "ancestry",
        "stratum",
        "freq_class",
        "cumulative_count",
    ]
    for key in buckets:
        buckets[key].sort()

    with open_text(out_path, "wt") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for k, sample in enumerate(order):
            for (sk, fc), ranks in buckets.items():
                w.writerow(
                    {
                        "sample_index": k,
                        "sample_id": sample,
                        "ancestry": ancestry.get(sample, "oth"),
                        "stratum": sk,
                        "freq_class": fc,
                        "cumulative_count": bisect.bisect_right(ranks, k),
                    }
                )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--carriers", required=True)
    p.add_argument("--sample-ancestry", required=True)
    p.add_argument(
        "--out",
        help="Single-stratum output TSV (required unless --out-prefix is set)",
    )
    p.add_argument(
        "--out-prefix",
        help="If set with --strata-list, write <prefix>.discovery[.region|.cadd].tsv",
    )
    p.add_argument(
        "--strata",
        default="none",
        choices=["none", "region", "cadd", "region_cadd"],
        help="Single stratification (default: none)",
    )
    p.add_argument(
        "--strata-list",
        default="",
        help="Comma-separated strata to emit in one pass (e.g. none,region,cadd)",
    )
    p.add_argument(
        "--include-sources",
        default="main",
        help="Comma-separated source_vcf values to include (default: main)",
    )
    args = p.parse_args()

    strata_list = [s.strip() for s in args.strata_list.split(",") if s.strip()]
    if not strata_list:
        strata_list = [args.strata]
    if args.out_prefix:
        out_map = {
            "none": Path(f"{args.out_prefix}.discovery.tsv"),
            "region": Path(f"{args.out_prefix}.discovery.region.tsv"),
            "cadd": Path(f"{args.out_prefix}.discovery.cadd.tsv"),
            "region_cadd": Path(f"{args.out_prefix}.discovery.region_cadd.tsv"),
        }
    else:
        if not args.out:
            raise SystemExit("Provide --out or --out-prefix")
        if len(strata_list) != 1:
            raise SystemExit("--out requires a single --strata value")
        out_map = {strata_list[0]: Path(args.out)}

    ancestry = load_ancestry(args.sample_ancestry)
    order = order_samples_by_ancestry(ancestry)
    include = set(args.include_sources.split(","))
    sample_rank = {s: i for i, s in enumerate(order)}

    buckets_by_strata: dict[str, dict[tuple[str, str], list[int]]] = {
        s: defaultdict(list) for s in strata_list
    }

    n = 0
    n_used = 0
    with open_text(args.carriers, "rt") as carr_fh:
        carr_header = carr_fh.readline().rstrip("\n").split("\t")
        i_id = carr_header.index("id")
        i_c = carr_header.index("carriers")
        sites_iter = iter_sites(args.sites)
        for site, carr_line in zip(sites_iter, carr_fh):
            n += 1
            if n % 500000 == 0:
                print(f"[ebert_discovery] {n:,} rows", file=sys.stderr, flush=True)
            parts = carr_line.rstrip("\n").split("\t")
            carr_id = parts[i_id] if i_id < len(parts) else ""
            if site["id"] != carr_id:
                raise RuntimeError(
                    f"sites/carriers order mismatch at row {n}: "
                    f"site={site['id']} carriers={carr_id}"
                )
            if site["source_vcf"] not in include:
                continue
            carriers_field = parts[i_c] if i_c < len(parts) else ""
            first = first_carrier_rank(carriers_field, sample_rank)
            if first is None:
                continue
            n_used += 1
            freq = site["freq_class"]
            region = site.get("region_class", "non_repetitive")
            cadd = site.get("cadd_sv_bin", "unscored")
            for strata in strata_list:
                sk = stratum_key(strata, region, cadd)
                buckets_by_strata[strata][(sk, freq)].append(first)

        extra_site = next(sites_iter, None)
        extra_carr = next(carr_fh, None)
        if extra_site is not None or extra_carr is not None:
            raise RuntimeError("sites and carriers TSV lengths differ")

    print(
        f"[ebert_discovery] scanned {n:,} rows; used {n_used:,} for discovery",
        file=sys.stderr,
        flush=True,
    )

    for strata in strata_list:
        out_path = out_map[strata]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_discovery(out_path, order, ancestry, buckets_by_strata[strata])
        print(f"[ebert_discovery] wrote {out_path}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
