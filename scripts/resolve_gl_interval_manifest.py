#!/usr/bin/env python3
"""Resolve DeepVariant joint-call VCF shards from GL_INTERVAL_set (TSV or Terra firecloud)."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd

DEFAULT_NAMESPACE = "allofus-drc-wgs-LR-prodData"
DEFAULT_WORKSPACE = "AoU_DRC_LongReads_PhaseTwo_Storage"
DEFAULT_ENTITY_TYPE = "GL_INTERVAL_set"
DEFAULT_ID_COLUMN = "entity:GL_INTERVAL_set_id"
DEFAULT_URI_COLUMN = "VCF"
DEFAULT_IDX_COLUMN = "VCF_idx"


def _entities_to_rows(entities: list[dict[str, Any]], *, id_column: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for ent in entities:
        entity_id = str(ent.get("name") or "").strip()
        if not entity_id:
            continue
        attrs = ent.get("attributes") or {}
        row = {id_column: entity_id}
        for key, value in attrs.items():
            row[str(key)] = "" if value is None else str(value).strip()
        rows.append(row)
    return rows


def load_gl_interval_manifest_tsv(path: Path, *, id_column: str = DEFAULT_ID_COLUMN) -> pd.DataFrame:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No rows in {path}")
    manifest = pd.DataFrame(rows, dtype=str)
    if id_column not in manifest.columns:
        candidates = [c for c in manifest.columns if c.startswith("entity:") or c.endswith("_id")]
        if not candidates:
            raise ValueError(f"No entity id column in {path}; columns={list(manifest.columns)}")
        manifest = manifest.rename(columns={candidates[0]: id_column})
    return manifest


def fetch_gl_interval_manifest_firecloud(
    namespace: str,
    workspace: str,
    entity_type: str = DEFAULT_ENTITY_TYPE,
    *,
    id_column: str = DEFAULT_ID_COLUMN,
) -> pd.DataFrame:
    from firecloud import api as fapi

    resp = fapi.get_entities(namespace, workspace, entity_type)
    if resp.status_code != 200:
        raise RuntimeError(
            f"firecloud get_entities failed ({resp.status_code}): {resp.text[:500]}"
        )
    entities = resp.json()
    if not isinstance(entities, list):
        raise RuntimeError(f"Unexpected firecloud payload type: {type(entities)}")
    rows = _entities_to_rows(entities, id_column=id_column)
    if not rows:
        raise RuntimeError(f"No entities returned for {namespace}/{workspace}/{entity_type}")
    return pd.DataFrame(rows, dtype=str)


def normalize_gl_interval_manifest(
    manifest: pd.DataFrame,
    *,
    id_column: str = DEFAULT_ID_COLUMN,
    uri_column: str = DEFAULT_URI_COLUMN,
    idx_column: str = DEFAULT_IDX_COLUMN,
    autosomes_only: bool = True,
    source_label: str = "manifest",
) -> pd.DataFrame:
    if id_column not in manifest.columns:
        candidates = [c for c in manifest.columns if c.startswith("entity:") or c.endswith("_id")]
        if not candidates:
            raise ValueError(f"No entity id column; columns={list(manifest.columns)}")
        id_column = candidates[0]
    if uri_column not in manifest.columns:
        raise ValueError(f"Missing {uri_column}; columns={list(manifest.columns)}")

    out = manifest.rename(columns={id_column: "interval_id"}).copy()
    out["interval_id"] = out["interval_id"].astype(str).str.strip()
    out[uri_column] = out[uri_column].astype(str).str.strip()
    out = out.replace({"": pd.NA, "nan": pd.NA})
    out = out.dropna(subset=[uri_column])

    if idx_column not in out.columns:
        out[idx_column] = out[uri_column].astype(str) + ".tbi"
    else:
        out[idx_column] = out[idx_column].astype(str).str.strip()
        missing_idx = out[idx_column].isna() | (out[idx_column] == "")
        out.loc[missing_idx, idx_column] = out.loc[missing_idx, uri_column].astype(str) + ".tbi"

    if autosomes_only:
        autosomes = {f"chr{i}" for i in range(1, 23)}
        out = out.loc[out["interval_id"].isin(autosomes)].copy()

    if out.empty:
        raise ValueError(f"No autosomal VCF rows in {source_label}")
    if not out["interval_id"].is_unique:
        dupes = out.loc[out["interval_id"].duplicated(), "interval_id"].tolist()
        raise ValueError(f"Duplicate interval ids in {source_label}: {dupes[:5]}")
    return out.reset_index(drop=True)


def resolve_gl_interval_manifest(
    *,
    namespace: str = DEFAULT_NAMESPACE,
    workspace: str = DEFAULT_WORKSPACE,
    entity_type: str = DEFAULT_ENTITY_TYPE,
    manifest_tsv: Path | None = None,
    entities_json: Path | None = None,
    from_firecloud: bool = False,
    id_column: str = DEFAULT_ID_COLUMN,
    uri_column: str = DEFAULT_URI_COLUMN,
    idx_column: str = DEFAULT_IDX_COLUMN,
    autosomes_only: bool = True,
) -> pd.DataFrame:
    if manifest_tsv is not None:
        manifest = load_gl_interval_manifest_tsv(manifest_tsv, id_column=id_column)
        source_label = str(manifest_tsv)
    elif entities_json is not None:
        entities = json.loads(entities_json.read_text())
        if not isinstance(entities, list):
            raise ValueError(f"Expected list in {entities_json}, got {type(entities)}")
        manifest = pd.DataFrame(_entities_to_rows(entities, id_column=id_column), dtype=str)
        source_label = str(entities_json)
    elif from_firecloud:
        manifest = fetch_gl_interval_manifest_firecloud(
            namespace, workspace, entity_type, id_column=id_column
        )
        source_label = f"{namespace}/{workspace}/{entity_type}"
    else:
        raise ValueError("Provide manifest_tsv, entities_json, or from_firecloud=True")

    return normalize_gl_interval_manifest(
        manifest,
        id_column=id_column,
        uri_column=uri_column,
        idx_column=idx_column,
        autosomes_only=autosomes_only,
        source_label=source_label,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--manifest-tsv", type=Path, help="Local TSV export of GL_INTERVAL_set")
    src.add_argument(
        "--from-firecloud",
        action="store_true",
        help="Fetch GL_INTERVAL_set via firecloud API",
    )
    src.add_argument("--entities-json", type=Path, help="Saved firecloud get_entities JSON")
    parser.add_argument("--namespace", default=DEFAULT_NAMESPACE)
    parser.add_argument("--workspace", default=DEFAULT_WORKSPACE)
    parser.add_argument("--entity-type", default=DEFAULT_ENTITY_TYPE)
    parser.add_argument("--out-tsv", type=Path)
    args = parser.parse_args()

    manifest = resolve_gl_interval_manifest(
        namespace=args.namespace,
        workspace=args.workspace,
        entity_type=args.entity_type,
        manifest_tsv=args.manifest_tsv,
        entities_json=args.entities_json,
        from_firecloud=args.from_firecloud,
    )
    print(manifest.to_csv(sep="\t", index=False))
    if args.out_tsv:
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(args.out_tsv, sep="\t", index=False)
        print(f"wrote {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
