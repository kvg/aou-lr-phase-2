#!/usr/bin/env python3
"""Sum 1 bp signed-length × region_class count TSVs and JSON summaries."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

BINS_HEADER = "chrom\tpartition\tsigned_len\tregion_class\tn_sites\n"

# Table 2 (tab:callset) cross-check, from the merged summary counts.
TABLE2_KEYS = {
    "snv_records": "SNVs (bcftools 'number of SNPs')",
    "ins_lt_small_max": "INS < 20 bp (GLnexus)",
    "del_lt_small_max": "DEL < 20 bp (GLnexus)",
    "main_ge20_DEL": "DEL >= 20 bp (main)",
    "main_ge50_DEL": "DEL >= 50 bp (main)",
    "main_ge20_INS": "INS >= 20 bp (main)",
    "main_ge50_INS": "INS >= 50 bp (main)",
    "large_records": "Large events > 10 kb",
    "bnd_records": "BND",
}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bins", nargs="+", required=True, help="Per-chrom and companion count TSVs")
    p.add_argument("--summaries", nargs="+", required=True, help="Matching summary JSONs")
    p.add_argument("--out-bins", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args()

    counts: Counter = Counter()
    for path in args.bins:
        with open(path, encoding="utf-8") as fh:
            header = fh.readline()
            if header != BINS_HEADER:
                raise SystemExit(
                    f"unexpected header in {path}: {header!r}\n"
                    f"Expected {BINS_HEADER!r} (summarize_variant_length_spectrum.py)."
                )
            for line in fh:
                _chrom, partition, signed_len, region, n_sites = line.rstrip("\n").split("\t")
                counts[(partition, int(signed_len), region)] += int(n_sites)

    out_bins = Path(args.out_bins)
    out_bins.parent.mkdir(parents=True, exist_ok=True)
    with out_bins.open("wt", encoding="utf-8") as out:
        out.write(BINS_HEADER)
        for (partition, slen, region), n in sorted(counts.items()):
            out.write(f"all\t{partition}\t{slen}\t{region}\t{n}\n")

    totals: Counter = Counter()
    chroms = []
    params: dict = {}
    for path in args.summaries:
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        chroms.append(d.get("chrom"))
        for k, v in (d.get("counts") or {}).items():
            totals[k] += int(v)
        for k in ("small_max_bp", "sv_min_bp", "apply_filters", "region_rule"):
            if k in d:
                params.setdefault(k, d[k])

    partition_totals: Counter = Counter()
    for (partition, _slen, _region), n in counts.items():
        partition_totals[partition] += n

    summary = {
        "chroms": chroms,
        "binning": "1bp_signed_len",
        **params,
        "n_length_rows": len(counts),
        "partition_totals": dict(sorted(partition_totals.items())),
        "table2_check": {label: int(totals.get(k, 0)) for k, label in TABLE2_KEYS.items()},
        "counts": {k: int(v) for k, v in sorted(totals.items())},
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
