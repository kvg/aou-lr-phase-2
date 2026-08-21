#!/usr/bin/env python3
"""
Build Ebert-style cumulative discovery tables.

Outputs a long TSV with cumulative counts after each ancestry-ordered sample,
optionally stratified by region_class and/or cadd_sv_bin.

Joins sites to carriers by variant ``id`` (order-independent). Streams sites
into discovery buckets so a multi-million-site table does not sit in RAM as
row dicts; only a compact id→first-carrier-rank map is retained.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import open_text, order_samples_by_ancestry  # noqa: E402


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


def load_first_carrier_rank(path: str, sample_rank: dict[str, int]) -> dict[str, int]:
    """Map site id -> earliest ancestry-ordered carrier rank (ints only)."""
    out: dict[str, int] = {}
    n = 0
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index("id")
        i_c = header.index("carriers")
        for line in fh:
            n += 1
            if n % 500000 == 0:
                print(f"[ebert_discovery] carriers {n:,}", file=sys.stderr, flush=True)
            parts = line.rstrip("\n").split("\t")
            if i_id >= len(parts):
                continue
            carriers_field = parts[i_c] if i_c < len(parts) else ""
            first = first_carrier_rank(carriers_field, sample_rank)
            if first is not None:
                out[parts[i_id]] = first
    print(f"[ebert_discovery] first-rank map: {len(out):,} / {n:,}", file=sys.stderr, flush=True)
    return out


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
    first_rank_by_id = load_first_carrier_rank(args.carriers, sample_rank)

    buckets_by_strata: dict[str, dict[tuple[str, str], list[int]]] = {
        s: defaultdict(list) for s in strata_list
    }

    n = 0
    n_used = 0
    with open_text(args.sites, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for required in ("id", "source_vcf", "freq_class"):
            if required not in col:
                raise SystemExit(f"sites TSV missing column: {required}")
        i_id = col["id"]
        i_src = col["source_vcf"]
        i_freq = col["freq_class"]
        i_region = col.get("region_class")
        i_cadd = col.get("cadd_sv_bin")
        for line in fh:
            n += 1
            if n % 500000 == 0:
                print(f"[ebert_discovery] sites {n:,}", file=sys.stderr, flush=True)
            parts = line.rstrip("\n").split("\t")
            if parts[i_src] not in include:
                continue
            first = first_rank_by_id.get(parts[i_id])
            if first is None:
                continue
            n_used += 1
            freq = parts[i_freq]
            region = parts[i_region] if i_region is not None and i_region < len(parts) else "non_repetitive"
            cadd = parts[i_cadd] if i_cadd is not None and i_cadd < len(parts) else "unscored"
            for strata in strata_list:
                sk = stratum_key(strata, region, cadd)
                buckets_by_strata[strata][(sk, freq)].append(first)

    # Free the join map before writing long discovery TSVs.
    first_rank_by_id.clear()

    print(
        f"[ebert_discovery] scanned {n:,} sites; used {n_used:,} for discovery",
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
