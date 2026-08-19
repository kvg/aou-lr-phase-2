#!/usr/bin/env python3
"""
Attach CADD-SV scores to the site table.

Accepts a CADD-SV v2 score TSV. Scores are joined by the normalized
coordinate key (chrom, 0-based start, 1-based end, SVTYPE); an ID join is
also supported for legacy score tables. Sites without scores stay unscored.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import (  # noqa: E402
    cadd_sv_bin,
    iter_sites,
    open_text,
    site_row,
    write_site_header,
)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--scores", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--id-col", default="id")
    p.add_argument("--phred-col", default="CADD-SV_PHRED")
    args = p.parse_args()

    scores_by_id: dict[str, float] = {}
    scores_by_coordinate: dict[tuple[str, int, int, str], float] = {}
    with open_text(args.scores, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # tolerate alternate headers
        id_col = args.id_col
        phred_col = args.phred_col
        if id_col not in header:
            for cand in ("id", "ID", "name", "Name"):
                if cand in header:
                    id_col = cand
                    break
        if phred_col not in header:
            for cand in ("CADD-SV_PHRED", "CADD_SV_PHRED", "PHRED", "phred"):
                if cand in header:
                    phred_col = cand
                    break
        i_ph = header.index(phred_col)
        i_id = header.index(id_col) if id_col in header else None

        def locate(*candidates: str):
            lowered = {name.lstrip("#").lower(): i for i, name in enumerate(header)}
            for candidate in candidates:
                if candidate.lower() in lowered:
                    return lowered[candidate.lower()]
            return None

        i_chrom = locate("chrom", "chr", "chromosome")
        i_start = locate("start", "pos", "position")
        i_end = locate("end")
        i_type = locate("type", "svtype")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                phred = float(parts[i_ph])
            except (ValueError, IndexError):
                continue
            if i_id is not None and i_id < len(parts):
                scores_by_id[parts[i_id]] = phred
            if None not in (i_chrom, i_start, i_end, i_type):
                try:
                    scores_by_coordinate[
                        (
                            parts[i_chrom],
                            int(parts[i_start]),
                            int(parts[i_end]),
                            parts[i_type].split(":")[0],
                        )
                    ] = phred
                except (ValueError, IndexError):
                    continue

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(args.out, "wt") as out:
        write_site_header(out)
        for s in iter_sites(args.sites):
            n += 1
            start0 = int(s["pos"]) - 1
            end1 = int(s["end"]) if s.get("end") not in (None, "") else start0 + abs(int(s["svlen"] or 0))
            ph = scores_by_coordinate.get((s["chrom"], start0, end1, s["svtype"]))
            if ph is None:
                ph = scores_by_id.get(s["id"])
            s["cadd_sv_phred"] = "" if ph is None else ph
            s["cadd_sv_bin"] = cadd_sv_bin(ph)
            out.write(site_row(s))
    print(f"[attach_caddsv] scored {n:,} sites", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
