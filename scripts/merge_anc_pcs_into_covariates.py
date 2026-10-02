#!/usr/bin/env python3
"""Add ancestry-export PCs to a covariates release.

Joins PC1–PC16 from ``anc_df.csv.gz`` as ``anc_PC1``–``anc_PC16``. Those are the
coordinates that also appear in the circulated long-read metadata
(``merged_all_df.csv.gz``). They do not replace full-training ``PC1``–``PC32``
or long-read ``lr_PC*``.

Writes a new covariates file and data dictionary. Does not modify the inputs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

ANC_PC_COLS = [f"anc_PC{i}" for i in range(1, 17)]
SRC_PC_COLS = [f"PC{i}" for i in range(1, 17)]


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def load_anc_pcs(path: Path, wanted: set[str]) -> dict[str, tuple[str, ...]]:
    """Return person_id -> 16 PC strings for IDs in ``wanted``.

    Blank coordinates are stored as empty strings. Duplicate IDs must agree.
    """
    found: dict[str, tuple[str, ...]] = {}
    with _open_text(path) as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in ["person_id", *SRC_PC_COLS] if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path} missing columns: {', '.join(missing)}")
        for row in reader:
            person_id = row["person_id"].strip()
            if person_id not in wanted:
                continue
            values = tuple(row[c].strip() for c in SRC_PC_COLS)
            previous = found.get(person_id)
            if previous is not None and previous != values:
                raise SystemExit(f"{path} has conflicting PC rows for {person_id}")
            found[person_id] = values
    return found


def insert_columns(fieldnames: list[str]) -> list[str]:
    if "anc_PC1" in fieldnames or "has_anc_pcs" in fieldnames:
        raise SystemExit("covariates already contain anc_PC* / has_anc_pcs")
    if "PC32" not in fieldnames:
        raise SystemExit("covariates have no PC32 column to insert anc_PC* after")
    if "has_global_pcs" not in fieldnames:
        raise SystemExit("covariates have no has_global_pcs column")
    cols = list(fieldnames)
    at = cols.index("PC32") + 1
    cols[at:at] = ANC_PC_COLS
    at = cols.index("has_global_pcs") + 1
    cols.insert(at, "has_anc_pcs")
    return cols


def update_dictionary(src: Path, dest: Path, final_columns: list[str]) -> None:
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_col = {row["column"]: row for row in rows}
    for i, column in enumerate(ANC_PC_COLS, start=1):
        by_col[column] = {
            "column": column,
            "source": "anc_df.csv.gz PC"
            + str(i)
            + " (same values in merged_all_df.csv.gz)",
            "notes": (
                f"ancestry-export PC{i} from the circulated long-read metadata; "
                f"does not replace full-training PC{i} or long-read lr_PC{i}"
            ),
        }
    by_col["has_anc_pcs"] = {
        "column": "has_anc_pcs",
        "source": "derived",
        "notes": "True if anc_PC1 is non-missing after the anc_df join",
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
            writer.writerow(
                by_col.get(column, {"column": column, "source": "", "notes": ""})
            )


def merge(covariates: Path, anc_csv: Path, out_csv: Path, dictionary: Path, out_dictionary: Path) -> None:
    with _open_text(covariates) as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "research_id" not in reader.fieldnames:
            raise SystemExit(f"{covariates} has no research_id column")
        base_fields = list(reader.fieldnames)
        cov_rows = list(reader)

    ids = [row["research_id"] for row in cov_rows]
    if len(ids) != len(set(ids)):
        raise SystemExit(f"{covariates} has duplicate research_id values")

    anc = load_anc_pcs(anc_csv, set(ids))
    out_fields = insert_columns(base_fields)
    n_hit = 0
    n_partial = 0
    for row in cov_rows:
        values = anc.get(row["research_id"])
        if values is None or values[0] == "":
            for column in ANC_PC_COLS:
                row[column] = ""
            row["has_anc_pcs"] = "False"
            if values is not None and any(values):
                n_partial += 1
            continue
        if any(value == "" for value in values):
            n_partial += 1
        for column, value in zip(ANC_PC_COLS, values):
            row[column] = value
        row["has_anc_pcs"] = "True"
        n_hit += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(out_csv).endswith(".gz") else open
    with opener(out_csv, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cov_rows)

    update_dictionary(dictionary, out_dictionary, out_fields)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_dictionary}")
    print(f"rows: {len(cov_rows):,}")
    print(f"has_anc_pcs: {n_hit:,}")
    print(f"partial anc PC rows (not flagged): {n_partial:,}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariates",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.csv.gz"),
    )
    parser.add_argument(
        "--anc-csv",
        type=Path,
        default=Path("tractor_mix/resources/legacy_covariates/anc_df.csv.gz"),
    )
    parser.add_argument(
        "--data-dictionary",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.data_dictionary.tsv"),
    )
    parser.add_argument("--out-covariates", type=Path, default=Path("covariates.v7.csv.gz"))
    parser.add_argument(
        "--out-data-dictionary",
        type=Path,
        default=Path("covariates.v7.data_dictionary.tsv"),
    )
    args = parser.parse_args()
    merge(
        args.covariates,
        args.anc_csv,
        args.out_covariates,
        args.data_dictionary,
        args.out_data_dictionary,
    )


if __name__ == "__main__":
    main()
