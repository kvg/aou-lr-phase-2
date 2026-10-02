"""Cohort CSV + v9 CRAM manifest helpers."""

from __future__ import annotations

import csv
import io
import os
import random
from pathlib import Path

BATCH_FIELDS = ("sample_id", "cram", "crai", "ref_fa", "ref_fai", "catalog", "sex")
LOCITYPER_FIELDS = (
    "sample_id",
    "cram",
    "crai",
    "ref_fa",
    "ref_fai",
    "counts_jf",
    "bed",
    "db_tar",
    "sex",
)


def find_v9_manifest() -> Path | None:
    env = (
        os.environ.get("EH_CRAM_MANIFEST")
        or os.environ.get("LOCITYPER_CRAM_MANIFEST")
        or ""
    ).strip()
    if env:
        p = Path(env)
        if p.is_file():
            return p.resolve()
    names = [
        Path("workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv"),
        Path("v9/wgs/cram/manifest.csv"),
    ]
    here = Path.cwd().resolve()
    roots = [here, *here.parents, Path("/home/jupyter"), Path("/home/jupyter/workspace")]
    for root in roots:
        for rel in names:
            cand = root / rel
            if cand.is_file():
                return cand.resolve()
    ws = Path("workspace")
    if ws.is_dir():
        hits = list(ws.glob("**/v9/wgs/cram/manifest.csv"))
        if hits:
            return hits[0].resolve()
    return None


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [{k: (v or "").strip() for k, v in rec.items()} for rec in csv.DictReader(fh)]


def write_csv(path: Path, rows: list[dict[str, str]], *, fields: tuple[str, ...] = BATCH_FIELDS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_tasks_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("no task rows")
    fields = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def records_from_manifest(
    path: Path,
    *,
    ref_fa: str,
    ref_fai: str,
    catalog: str,
    sex: str = "female",
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            sid = (rec.get("person_id") or rec.get("sample_id") or "").strip()
            cram = (rec.get("cram_uri") or rec.get("cram") or "").strip()
            crai = (rec.get("cram_index_uri") or rec.get("crai") or "").strip()
            if not sid or not cram or not crai:
                continue
            rows.append(
                {
                    "sample_id": sid,
                    "cram": cram,
                    "crai": crai,
                    "ref_fa": ref_fa,
                    "ref_fai": ref_fai,
                    "catalog": catalog,
                    "sex": (rec.get("sex") or sex).strip() or sex,
                }
            )
    return rows


def sample_pilot(rows: list[dict[str, str]], n: int, seed: int) -> list[dict[str, str]]:
    if n <= 0:
        raise ValueError("pilot n must be > 0")
    if n >= len(rows):
        return list(rows)
    rng = random.Random(seed)
    return rng.sample(list(rows), n)


def shards(rows: list[dict[str, str]], size: int) -> list[tuple[int, list[dict[str, str]]]]:
    if size <= 0:
        raise ValueError("shard size must be > 0")
    out: list[tuple[int, list[dict[str, str]]]] = []
    for i in range(0, len(rows), size):
        out.append((len(out), rows[i : i + size]))
    return out


SAMPLE_ID_COLUMNS = ("sample_id", "person_id", "research_id", "participant_id")
CRAM_COLUMNS = ("cram", "cram_uri", "cram_path")
CRAI_COLUMNS = ("crai", "cram_index_uri", "crai_uri", "crai_path", "cram_index")
SEX_COLUMNS = ("sex", "sex_at_birth", "inferred_sex", "gender")


def _pick(header: list[str], names: tuple[str, ...]) -> str | None:
    lower = {h.strip().lower(): h for h in header}
    for name in names:
        if name in lower:
            return lower[name]
    return None


def normalize_sex(raw: str) -> tuple[str, bool]:
    """ExpansionHunter --sex value and whether raw was recognised. Unknown values become female."""
    key = (raw or "").strip().lower().replace(" ", "").replace("_", "")
    if key in {"male", "m", "1", "xy"}:
        return "male", True
    if key in {"female", "f", "2", "xx"}:
        return "female", True
    return "female", False


def read_sample_table(text: str) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Parse a sample_id / CRAM / CRAI / sex table (CSV or TSV, flexible column names).

    Returns (records, report). Records carry sample_id, cram, crai, sex (male or
    female) and sex_raw. The report counts raw sex values, rows dropped and why,
    and duplicate sample ids. Nothing is silently fixed except unknown sex.
    """
    first = text.splitlines()[0] if text else ""
    delimiter = "\t" if first.count("\t") > first.count(",") else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    header = list(reader.fieldnames or [])
    cols = {
        "sample_id": _pick(header, SAMPLE_ID_COLUMNS),
        "cram": _pick(header, CRAM_COLUMNS),
        "crai": _pick(header, CRAI_COLUMNS),
        "sex": _pick(header, SEX_COLUMNS),
    }
    missing = [k for k, v in cols.items() if v is None]
    if missing:
        raise ValueError(f"sample table is missing columns for {missing}; header is {header}")
    records: list[dict[str, str]] = []
    sex_raw: dict[str, int] = {}
    dropped: dict[str, int] = {}
    seen: set[str] = set()
    dups: list[str] = []
    for row in reader:
        sid = (row.get(cols["sample_id"]) or "").strip()
        cram = (row.get(cols["cram"]) or "").strip()
        crai = (row.get(cols["crai"]) or "").strip()
        raw = (row.get(cols["sex"]) or "").strip()
        if not sid:
            dropped["no sample id"] = dropped.get("no sample id", 0) + 1
            continue
        if not cram.startswith("gs://") or not crai.startswith("gs://"):
            dropped["CRAM or CRAI not gs://"] = dropped.get("CRAM or CRAI not gs://", 0) + 1
            continue
        if sid in seen:
            dups.append(sid)
            continue
        seen.add(sid)
        sex, _ = normalize_sex(raw)
        sex_raw[raw or "(blank)"] = sex_raw.get(raw or "(blank)", 0) + 1
        records.append({"sample_id": sid, "cram": cram, "crai": crai, "sex": sex, "sex_raw": raw})
    unrecognized = {k: v for k, v in sex_raw.items() if not normalize_sex("" if k == "(blank)" else k)[1]}
    report = {
        "columns": cols,
        "n": len(records),
        "sex_raw": dict(sorted(sex_raw.items(), key=lambda kv: -kv[1])),
        "sex_unrecognized_as_female": unrecognized,
        "dropped": dropped,
        "duplicates": dups,
    }
    return records, report
