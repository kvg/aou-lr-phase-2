#!/usr/bin/env python3
"""Expand one FelixGenome run-table row into a batch of rows over a phenotype list.

A genome-wide scan over the full phecode set is 71P + 92 Cromwell calls, which
at P = 943 is 67,045 in a single workflow. Splitting the phenotype list across
submissions costs nothing — the GRM and all 22 chromosome packs call-cache
across batches — and buys workflows that can be watched and restarted.

Batches share one pheno_cov.tsv. Only `selected_phenotypes` differs between
rows: fit_felix_null.R reads just the named phenotype column, and the
summariser reads pheno_cov with usecols, so a table carrying every phecode
column costs a batch nothing.

Writes the entity and membership TSVs to upload to Terra, and the per-batch
phenotype files to copy to the bucket.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

CHROMS = [f"chr{i}" for i in range(1, 23)]
ID_COL = "entity:flare_lai_prod_set_id"
DEFAULT_SET_TYPE = "flare_lai_prod_set"
DEFAULT_MEMBER_TYPE = "flare_lai_prod"
DEFAULT_MEMBER_PREFIX = "aou_lr_phase2_v1"


def read_entity_tsv(path: Path) -> tuple[list[str], dict[str, dict], str]:
    """Columns, rows keyed by entity id, and the id column's own name.

    The id column is read from the header rather than assumed, so a table can be
    cloned into a differently named set type (--set-type) while still being read
    from whatever Terra exported.
    """
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = list(rdr.fieldnames or [])
        if not cols or not (cols[0].startswith("entity:") and cols[0].endswith("_id")):
            raise SystemExit(
                f"{path} is not an entity TSV: first column must look like "
                f"entity:<type>_id, got {cols[0] if cols else '(empty)'!r}")
        id_col = cols[0]
        rows = {r[id_col]: r for r in rdr}
    return cols, rows, id_col


def batches(names: list[str], size: int) -> list[list[str]]:
    """Split into the fewest batches of at most `size`, as evenly as possible.

    `size` is a ceiling, not a quota. Slicing at a fixed width leaves whatever
    is left over in the last batch, and 951 phenotypes at 95 leaves a batch of
    one -- a whole submission, with its 92 per-run calls, for a single
    phenotype. Spreading the remainder costs nothing and cannot produce that.
    """
    n_batches = math.ceil(len(names) / size)
    base, extra = divmod(len(names), n_batches)
    out, start = [], 0
    for i in range(n_batches):
        width = base + (1 if i < extra else 0)
        out.append(names[start:start + width])
        start += width
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phenotypes", type=Path, required=True,
                    help="selected_phenotypes.txt for the whole run, one name per line. "
                         "Order is the contract: null arrays are indexed against it.")
    ap.add_argument("--template", type=Path, required=True,
                    help="entity TSV holding the row to clone")
    ap.add_argument("--template-row", required=True, help="id of the row to clone")
    ap.add_argument("--row-prefix", default=None,
                    help="prefix for the generated row ids [default: the template row id]")
    ap.add_argument("--batch-size", type=int, default=95,
                    help="phenotypes per submission [default: %(default)s, "
                         "about 7,000 calls per workflow]")
    ap.add_argument("--bucket", required=True,
                    help="gs://BUCKET the per-batch phenotype files will be copied to")
    ap.add_argument("--phenotype-prefix", default="felix/phenotypes",
                    help="path under the bucket for the phenotype files [default: %(default)s]")
    ap.add_argument("--pheno-cov", default=None,
                    help="override the template's pheno_cov cell (the new table with every "
                         "phecode column); shorthand for --set pheno_cov=...")
    ap.add_argument("--set", dest="overrides", action="append", default=[],
                    metavar="COLUMN=VALUE",
                    help="override any column on every generated row, repeatable. The "
                         "column must already exist in the template.")
    ap.add_argument("--name-rows-by-phenotype", action="store_true",
                    help="name each row <prefix>.<phenotype> instead of <prefix>.batchNNN. "
                         "Only valid with --batch-size 1; the point of one row per "
                         "phenotype is being able to find the row for a failing "
                         "phenotype, which batch437 does not give you.")
    ap.add_argument("--set-type", default=DEFAULT_SET_TYPE,
                    help="Terra entity type to WRITE (default %(default)s). Use a new name "
                         "to clone an existing table into a separate one; the method config "
                         "then needs only its root entity type changed, because the column "
                         "names and the ${this.<member>s.*} references are unchanged.")
    ap.add_argument("--member-type", default=DEFAULT_MEMBER_TYPE,
                    help="membership column / member entity type (default %(default)s). "
                         "This is what ${this.<member>s.*} in the inputs JSON resolves "
                         "against, so changing it requires editing that JSON too.")
    ap.add_argument("--member-prefix", default=DEFAULT_MEMBER_PREFIX,
                    help="member ids are <prefix>.<chrom> (default %(default)s). Must match "
                         "the ids already in the table or the run gets different VCFs and "
                         "every cached call is lost.")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    names = [ln.strip() for ln in args.phenotypes.read_text().splitlines() if ln.strip()]
    if not names:
        raise SystemExit(f"{args.phenotypes} is empty")
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise SystemExit(f"duplicate phenotype names: {sorted(dupes)[:5]}")
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be at least 1")

    cols, rows, src_id_col = read_entity_tsv(args.template)
    out_id_col = f"entity:{args.set_type}_id"
    if args.template_row not in rows:
        raise SystemExit(f"row {args.template_row!r} not in {args.template}; "
                         f"have {sorted(rows)}")
    template = rows[args.template_row]
    if not any((template[c] or "").strip() for c in cols if c != src_id_col):
        raise SystemExit(f"row {args.template_row!r} has no values — it is probably one of "
                         "the rows Terra auto-creates from a submission")

    overrides: dict[str, str] = {}
    for item in args.overrides:
        if "=" not in item:
            raise SystemExit(f"--set expects COLUMN=VALUE, got {item!r}")
        key, value = item.split("=", 1)
        key = key.strip()
        if key not in cols:
            raise SystemExit(f"--set {key}: not a column in {args.template}; have {cols}")
        if key == src_id_col:
            raise SystemExit(f"--set cannot change {src_id_col}; use --row-prefix")
        overrides[key] = value
    if args.pheno_cov:
        overrides.setdefault("pheno_cov", args.pheno_cov)

    if args.name_rows_by_phenotype and args.batch_size != 1:
        raise SystemExit("--name-rows-by-phenotype requires --batch-size 1")
    prefix = args.row_prefix or args.template_row
    bucket = args.bucket.rstrip("/")
    pheno_dir = args.phenotype_prefix.strip("/")
    groups = batches(names, args.batch_size)
    width = max(2, len(str(len(groups))))

    out = args.out_dir
    (out / "phenotypes").mkdir(parents=True, exist_ok=True)

    new_rows = []
    for i, group in enumerate(groups, start=1):
        rid = (f"{prefix}.{group[0]}" if args.name_rows_by_phenotype
               else f"{prefix}.batch{i:0{width}d}")
        fname = f"{rid}.txt"
        (out / "phenotypes" / fname).write_text("\n".join(group) + "\n")
        row = dict(template)
        row.pop(src_id_col, None)
        row[out_id_col] = rid
        row["selected_phenotypes"] = f"{bucket}/{pheno_dir}/{fname}"
        row.update(overrides)
        if "run_description" in row:
            row["run_description"] = (
                f"Phenotype {group[0]} ({i} of {len(groups)}) from "
                f"{args.template_row}: {71 * len(group) + 92:,} calls."
                if args.name_rows_by_phenotype else
                f"Batch {i} of {len(groups)} from {args.template_row}: "
                f"{len(group)} phenotypes ({group[0]}..{group[-1]}), "
                f"{71 * len(group) + 92:,} calls.")
        new_rows.append(row)

    ent = out / "felix_runs.entity.tsv"
    with ent.open("w", newline="") as fh:
        out_cols = [out_id_col] + [c for c in cols if c != src_id_col]
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(new_rows)

    mem = out / "felix_runs.membership.tsv"
    with mem.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"membership:{args.set_type}_id", args.member_type])
        for row in new_rows:
            for c in CHROMS:
                w.writerow([row[out_id_col], f"{args.member_prefix}.{c}"])

    written = {out / "phenotypes" / f"{r[out_id_col]}.txt" for r in new_rows}
    stale = sorted(p for p in (out / "phenotypes").glob("*.txt") if p not in written)
    if stale:
        print(f"note: {len(stale)} file(s) in {out / 'phenotypes'} are from an earlier "
              f"run and are not part of this table: "
              + ", ".join(p.name for p in stale[:3]) + (" ..." if len(stale) > 3 else ""))

    total = sum(71 * len(g) + 92 for g in groups)
    print(f"{len(names)} phenotypes -> {len(groups)} batches of up to {args.batch_size}")
    print(f"  {ent}")
    print(f"  {mem}  ({len(new_rows) * len(CHROMS)} rows)")
    print(f"  calls: {71 * len(names) + 92:,} as one workflow, "
          f"{total:,} across {len(groups)} ({math.ceil(total / len(groups)):,} each)")
    print(f"\ngsutil -m cp {out / 'phenotypes'}/*.txt {bucket}/{pheno_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
