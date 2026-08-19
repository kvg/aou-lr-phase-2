#!/usr/bin/env python3
"""
Build Ebert-style cumulative discovery tables.

Outputs a long TSV with cumulative counts after each ancestry-ordered sample,
optionally stratified by region_class and/or cadd_sv_bin.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text, order_samples_by_ancestry  # noqa: E402


FREQ_ORDER = ["shared", "major", "polymorphic", "singleton"]


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


def load_first_carrier_rank(path: str, sample_rank: dict[str, int]) -> dict[str, int]:
    """Map site id -> earliest ancestry-ordered carrier rank."""
    out: dict[str, int] = {}
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index("id")
        i_c = header.index("carriers")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            carriers = parts[i_c] if i_c < len(parts) else ""
            if not carriers:
                continue
            best: int | None = None
            for sample in carriers.split(","):
                rank = sample_rank.get(sample)
                if rank is None:
                    continue
                if best is None or rank < best:
                    best = rank
            if best is not None:
                out[parts[i_id]] = best
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--carriers", required=True)
    p.add_argument("--sample-ancestry", required=True)
    p.add_argument("--out", required=True)
    p.add_argument(
        "--strata",
        default="none",
        choices=["none", "region", "cadd", "region_cadd"],
        help="Additional stratification beyond frequency class",
    )
    p.add_argument(
        "--include-sources",
        default="main",
        help="Comma-separated source_vcf values to include (default: main)",
    )
    args = p.parse_args()

    ancestry = load_ancestry(args.sample_ancestry)
    order = order_samples_by_ancestry(ancestry)
    include = set(args.include_sources.split(","))
    sample_rank = {s: i for i, s in enumerate(order)}
    first_rank_by_id = load_first_carrier_rank(args.carriers, sample_rank)

    site_meta = []
    for s in iter_sites(args.sites):
        if s["source_vcf"] not in include:
            continue
        first = first_rank_by_id.get(s["id"])
        if first is None:
            continue
        site_meta.append(
            {
                "first_rank": first,
                "freq_class": s["freq_class"],
                "region_class": s.get("region_class", "non_repetitive"),
                "cadd_sv_bin": s.get("cadd_sv_bin", "unscored"),
            }
        )

    def stratum_key(m):
        if args.strata == "none":
            return ("all",)
        if args.strata == "region":
            return (m["region_class"],)
        if args.strata == "cadd":
            return (m["cadd_sv_bin"],)
        return (m["region_class"], m["cadd_sv_bin"])

    # cumulative counts: after including samples 0..k
    # For each k, count sites with first_rank <= k, bucketed by freq_class and stratum
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_index",
        "sample_id",
        "ancestry",
        "stratum",
        "freq_class",
        "cumulative_count",
    ]
    # Aggregate efficiently: group sites by (stratum, freq_class) and first_rank
    buckets: dict[tuple, list[int]] = defaultdict(list)
    for m in site_meta:
        sk = "|".join(stratum_key(m))
        buckets[(sk, m["freq_class"])].append(m["first_rank"])
    for k in buckets:
        buckets[k].sort()

    with open_text(args.out, "wt") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for k, sample in enumerate(order):
            for (sk, fc), ranks in buckets.items():
                # count how many ranks <= k
                # linear scan is fine for modest N; use bisect for larger
                import bisect

                c = bisect.bisect_right(ranks, k)
                w.writerow(
                    {
                        "sample_index": k,
                        "sample_id": sample,
                        "ancestry": ancestry.get(sample, "oth"),
                        "stratum": sk,
                        "freq_class": fc,
                        "cumulative_count": c,
                    }
                )


if __name__ == "__main__":
    main()
