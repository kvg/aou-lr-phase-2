#!/usr/bin/env python3
"""Add eQTL, sQTL, and pQTL analysis-membership flags to a covariates release.

Reads unique research IDs from ``eqtl_ids.txt``, ``sqtl_ids.txt``, and
``pqtl_ids.txt`` inside ``covariates.tar.gz``. Duplicate lines in those files
are collapsed. Sets ``in_eqtl``, ``in_sqtl``, and ``in_pqtl`` to True when
``research_id`` is in the corresponding list, and False otherwise.

Does not add rows and does not replace ``has_rna`` or ``has_proteomics``.
Every listed ID must already be a covariates row. Writes a new covariates
file and data dictionary. Does not modify the inputs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from pathlib import Path

MEMBERSHIP = (
    ("in_eqtl", "eqtl_ids.txt", "eQTL"),
    ("in_sqtl", "sqtl_ids.txt", "sQTL"),
    ("in_pqtl", "pqtl_ids.txt", "pQTL"),
)
INSERT_AFTER = "has_proteomics"
HEADER_TOKENS = {"research_id", "sample_id", "id", "person_id"}


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def load_id_list(tar_path: Path, filename: str) -> tuple[set[str], int]:
    """Return unique IDs and the number of non-empty lines in ``filename``."""
    suffix = f"/{filename}"
    with tarfile.open(tar_path) as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and (member.name == filename or member.name.endswith(suffix))
        ]
        if len(members) != 1:
            names = ", ".join(member.name for member in members) or "(none)"
            raise SystemExit(f"{tar_path} must contain exactly one {filename}; found {names}")
        raw = archive.extractfile(members[0])
        if raw is None:
            raise SystemExit(f"{tar_path} could not read {members[0].name}")
        lines = raw.read().decode().splitlines()

    ids: set[str] = set()
    n_lines = 0
    for line in lines:
        token = line.strip().split()
        if not token:
            continue
        research_id = token[0]
        if n_lines == 0 and research_id.lower() in HEADER_TOKENS:
            continue
        n_lines += 1
        ids.add(research_id)
    if not ids:
        raise SystemExit(f"{tar_path} {filename} has no research IDs")
    return ids, n_lines


def insert_columns(fieldnames: list[str]) -> list[str]:
    present = [column for column, _, _ in MEMBERSHIP if column in fieldnames]
    if present:
        raise SystemExit(f"covariates already contain {', '.join(present)}")
    if INSERT_AFTER not in fieldnames:
        raise SystemExit(f"covariates have no {INSERT_AFTER} column to insert QTL flags after")
    columns = list(fieldnames)
    at = columns.index(INSERT_AFTER) + 1
    columns[at:at] = [column for column, _, _ in MEMBERSHIP]
    return columns


def update_dictionary(src: Path, dest: Path, final_columns: list[str], tar_path: Path) -> None:
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_col = {row["column"]: row for row in rows}
    for column, filename, label in MEMBERSHIP:
        by_col[column] = {
            "column": column,
            "source": f"{tar_path.name} {filename}",
            "notes": (
                f"True if research_id is listed in the {label} analysis sample list "
                f"({filename}); duplicate lines are collapsed to unique IDs; "
                "False otherwise. Does not replace has_rna or has_proteomics."
            ),
        }
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["column", "source", "notes"],
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for column in final_columns:
            writer.writerow(by_col.get(column, {"column": column, "source": "", "notes": ""}))


def merge(
    covariates: Path,
    ids_tar: Path,
    out_csv: Path,
    dictionary: Path,
    out_dictionary: Path,
) -> None:
    with _open_text(covariates) as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "research_id" not in reader.fieldnames:
            raise SystemExit(f"{covariates} has no research_id column")
        base_fields = list(reader.fieldnames)
        cov_rows = list(reader)

    ids = [row["research_id"] for row in cov_rows]
    if len(ids) != len(set(ids)):
        raise SystemExit(f"{covariates} has duplicate research_id values")
    id_set = set(ids)

    lists: dict[str, set[str]] = {}
    for column, filename, label in MEMBERSHIP:
        listed, n_lines = load_id_list(ids_tar, filename)
        missing = sorted(listed - id_set)
        if missing:
            preview = ", ".join(missing[:10])
            raise SystemExit(
                f"{label} list has {len(missing)} research IDs absent from {covariates}: {preview}"
            )
        lists[column] = listed
        print(f"{column}: {len(listed):,} unique IDs from {n_lines:,} lines")

    out_fields = insert_columns(base_fields)
    counts = {column: 0 for column, _, _ in MEMBERSHIP}
    for row in cov_rows:
        research_id = row["research_id"]
        for column, _, _ in MEMBERSHIP:
            flag = research_id in lists[column]
            row[column] = "True" if flag else "False"
            counts[column] += int(flag)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(out_csv).endswith(".gz") else open
    with opener(out_csv, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=out_fields,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(cov_rows)

    update_dictionary(dictionary, out_dictionary, out_fields, ids_tar)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_dictionary}")
    print(f"rows: {len(cov_rows):,}")
    for column, _, _ in MEMBERSHIP:
        print(f"{column}: {counts[column]:,}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covariates", type=Path, default=Path("covariates.v7.csv.gz"))
    parser.add_argument(
        "--ids-tar",
        type=Path,
        default=Path("../aou-lr-phase-2-manuscript/covariates.tar.gz"),
    )
    parser.add_argument(
        "--data-dictionary",
        type=Path,
        default=Path("covariates.v7.data_dictionary.tsv"),
    )
    parser.add_argument("--out-covariates", type=Path, default=Path("covariates.v8.csv.gz"))
    parser.add_argument(
        "--out-data-dictionary",
        type=Path,
        default=Path("covariates.v8.data_dictionary.tsv"),
    )
    args = parser.parse_args()
    merge(
        args.covariates,
        args.ids_tar,
        args.out_covariates,
        args.data_dictionary,
        args.out_data_dictionary,
    )


if __name__ == "__main__":
    main()
