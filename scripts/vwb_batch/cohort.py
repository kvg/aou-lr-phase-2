"""Cohort CSV + v9 CRAM manifest helpers."""

from __future__ import annotations

import csv
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
