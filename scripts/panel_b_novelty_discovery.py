#!/usr/bin/env python3
"""Cumulative Phase 2 SVs split by prior catalog: HPRC/HGSVC3, then Phase 1 (Fig. 2B).

Sites: main-callset DEL / INS with 50 <= |SVLEN| <= 10,000 (Table 2 rows), skipping
``suppressed_main_duplicate``. Classes, first match wins:

  * ``in_controls``   carried by >= 1 HPRC/HGSVC3 reference control
  * ``in_phase1``     Truvari-matched to a Phase 1 site (``--phase1-matched``)
  * ``new_singleton`` not in Phase 1, one carrier among participants
  * ``new_recurrent`` not in Phase 1, >= 2 carriers among participants

Participants follow ``--order`` (Fig. 2B discovery order). Each site enters the
curve at its first participant carrier; sites carried only by controls never
enter. Output is cumulative counts every ``--step``
participants (plus the last), with no sample ids, so it can leave the Workbench.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import open_text  # noqa: E402

CLASSES = ("in_controls", "in_phase1", "new_singleton", "new_recurrent")
AMBIGUOUS = -2


def _truthy(v: str) -> bool:
    return str(v or "").strip().lower() in {"1", "true", "t", "yes", "y"}


def load_order(path: str) -> tuple[list[str], list[str]]:
    ids, blocks = [], []
    with open_text(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ids.append(r["research_id"].strip())
            blocks.append(r.get("block", ""))
    return ids, blocks


def load_controls(path: str) -> set[str]:
    out: set[str] = set()
    with open_text(path, "rt") as fh:
        first = fh.readline()
        if "research_id" not in first.split("\t")[0]:
            return {s.strip() for s in [first, *fh] if s.strip()}
        fh.seek(0)
        rows = csv.DictReader(fh, delimiter="\t")
        for r in rows:
            if "is_reference_control" in r and not _truthy(r["is_reference_control"]):
                continue
            out.add((r.get("research_id") or next(iter(r.values()))).strip())
    return out


def load_matched(path: str) -> set[tuple[str, str, str]]:
    out: set[tuple[str, str, str]] = set()
    with open_text(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] != "chrom":
                out.add((parts[0], parts[1], parts[2]))
    return out


def load_carriers(path: str, rank: dict[str, int], controls: set[str], totals: Counter) -> dict[str, tuple[int, int, bool]]:
    """id -> (first participant rank or -1, n participant carriers, any control carrier); AMBIGUOUS in slot 0."""
    out: dict[str, tuple[int, int, bool]] = {}
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id, i_c = header.index("id"), header.index("carriers")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            vid = parts[i_id]
            if vid in out:
                totals["carriers_duplicate_id"] += 1
                out[vid] = (AMBIGUOUS, 0, False)
                continue
            field = parts[i_c] if i_c < len(parts) else ""
            first, n = None, 0
            has_control = False
            for s in field.split(",") if field else ():
                if s in controls:
                    has_control = True
                    continue
                r = rank.get(s)
                if r is None:
                    continue
                n += 1
                first = r if first is None or r < first else first
            if first is not None or has_control:
                out[vid] = (-1 if first is None else first, n, has_control)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True, help="Phase 2 main sites.tsv (AnnotateSvCallset)")
    p.add_argument("--carriers", required=True, help="Phase 2 main carriers.tsv (id, carriers, n_carriers)")
    p.add_argument("--order", required=True, help="research_id<TAB>block in discovery order")
    p.add_argument("--controls", required=True, help="control_sample_metadata.tsv or one id per line")
    p.add_argument("--phase1-matched", required=True, help="chrom<TAB>pos<TAB>id of Phase 2 sites matched to Phase 1")
    p.add_argument("--min-len", type=int, default=50)
    p.add_argument("--max-len", type=int, default=10_000)
    p.add_argument("--step", type=int, default=25)
    p.add_argument("--out", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args()

    order, blocks = load_order(args.order)
    rank = {s: i for i, s in enumerate(order)}
    controls = load_controls(args.controls)
    matched = load_matched(args.phase1_matched)
    totals: Counter = Counter()
    carriers = load_carriers(args.carriers, rank, controls, totals)

    buckets: dict[tuple[str, str], list[int]] = defaultdict(list)
    with open_text(args.sites, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        i_dup = col.get("suppressed_main_duplicate")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[col["source_vcf"]] != "main":
                continue
            svtype = parts[col["svtype"]].upper()
            if svtype not in ("DEL", "INS"):
                continue
            try:
                L = abs(int(parts[col["svlen"]]))
            except ValueError:
                continue
            if not args.min_len <= L <= args.max_len:
                continue
            if i_dup is not None and _truthy(parts[i_dup]):
                totals["skip_suppressed_duplicate"] += 1
                continue
            totals[f"sites_{svtype}"] += 1
            hit = carriers.get(parts[col["id"]])
            if hit is None:
                totals["skip_no_participant_carrier"] += 1
                continue
            first, n, has_control = hit
            if first == AMBIGUOUS:
                totals["skip_ambiguous_id"] += 1
                continue
            if first < 0:
                totals[f"controls_only_{svtype}"] += 1
                continue
            key = (parts[col["chrom"]], parts[col["pos"]], parts[col["id"]])
            if has_control:
                cls = "in_controls"
            elif key in matched:
                cls = "in_phase1"
            else:
                cls = "new_singleton" if n == 1 else "new_recurrent"
            totals[f"{cls}_{svtype}"] += 1
            buckets[(svtype, cls)].append(first)

    for v in buckets.values():
        v.sort()
    n_part = len(order)
    grid = list(range(args.step - 1, n_part, args.step))
    if not grid or grid[-1] != n_part - 1:
        grid.append(n_part - 1)
    with open_text(args.out, "wt") as out:
        out.write("n_participants\tblock\tsvtype\tclass\tcumulative\n")
        for k in grid:
            for svtype in ("INS", "DEL"):
                for cls in CLASSES:
                    c = bisect.bisect_right(buckets.get((svtype, cls), []), k)
                    out.write(f"{k + 1}\t{blocks[k]}\t{svtype}\t{cls}\t{c}\n")

    block_edges = {}
    for i, b in enumerate(blocks):
        block_edges.setdefault(b, [i + 1, i + 1])[1] = i + 1
    summary = {
        "n_participants": n_part,
        "n_controls": len(controls),
        "n_phase1_matched_keys": len(matched),
        "length_range": [args.min_len, args.max_len],
        "step": args.step,
        "block_ranges": block_edges,
        "counts": dict(sorted(totals.items())),
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary["counts"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
