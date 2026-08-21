#!/usr/bin/env python3
"""Resolve FLARE anc.vcf.gz URIs from aou_lr_chrom (TSV export or Terra firecloud table)."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_NAMESPACE = "allofus-drc-wgs-LR-prodData"
DEFAULT_WORKSPACE = "AoU_DRC_LongReads_PhaseTwo_Storage"
DEFAULT_ENTITY_TYPE = "aou_lr_chrom"
DEFAULT_URI_COLUMN = "model_chr_anc_vcf"


def chrom_from_entity_id(entity_id: str) -> str:
    # aou_lr_phase2_v1.chr22 -> chr22
    return entity_id.split(".")[-1] if entity_id else ""


def load_chrom_table_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            rid = row.get("entity:aou_lr_chrom_id") or row.get("aou_lr_chrom_id") or ""
            chrom = chrom_from_entity_id(rid)
            if chrom:
                rows[chrom] = {k: (v if v is not None else "") for k, v in row.items()}
        return rows


def load_chrom_table_entities(entities: list[dict[str, Any]]) -> dict[str, dict[str, str]]:
    """Convert firecloud get_entities JSON into chrom -> attribute map."""
    rows: dict[str, dict[str, str]] = {}
    for ent in entities:
        rid = str(ent.get("name") or "")
        chrom = chrom_from_entity_id(rid)
        if not chrom:
            continue
        attrs = ent.get("attributes") or {}
        flat = {str(k): "" if v is None else str(v) for k, v in attrs.items()}
        flat["aou_lr_chrom_id"] = rid
        rows[chrom] = flat
    return rows


def fetch_chrom_table_firecloud(
    namespace: str,
    workspace: str,
    entity_type: str = DEFAULT_ENTITY_TYPE,
) -> dict[str, dict[str, str]]:
    from firecloud import api as fapi

    resp = fapi.get_entities(namespace, workspace, entity_type)
    if resp.status_code != 200:
        raise SystemExit(
            f"firecloud get_entities failed ({resp.status_code}): {resp.text[:500]}"
        )
    entities = resp.json()
    if not isinstance(entities, list):
        raise SystemExit(f"Unexpected firecloud payload type: {type(entities)}")
    return load_chrom_table_entities(entities)


DEFAULT_AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def _uri_for(table: dict[str, dict[str, str]], chrom: str, uri_column: str) -> str:
    if chrom not in table:
        raise SystemExit(f"Chromosome not in table: {chrom}; have {sorted(table)}")
    u = (table[chrom].get(uri_column) or "").strip()
    if not u:
        raise SystemExit(
            f"{chrom}: empty {uri_column}. "
            f"global_anc={table[chrom].get('global_anc')!r}"
        )
    return u


def resolve_uris(
    table: dict[str, dict[str, str]],
    *,
    scan_chrom: str,
    grm_chroms: list[str],
    uri_column: str = DEFAULT_URI_COLUMN,
) -> dict[str, Any]:
    """Pilot helper: one scan VCF + GRM VCF list."""
    missing = [c for c in [scan_chrom, *grm_chroms] if c not in table]
    if missing:
        raise SystemExit(f"Chromosomes not in table: {missing}; have {sorted(table)}")

    return {
        "flare_vcf": _uri_for(table, scan_chrom, uri_column),
        "grm_vcfs": [_uri_for(table, c, uri_column) for c in grm_chroms],
        "scan_chrom": scan_chrom,
        "grm_chroms": grm_chroms,
        "uri_column": uri_column,
    }


def resolve_genome_uris(
    table: dict[str, dict[str, str]],
    *,
    scan_chroms: list[str],
    grm_chroms: list[str],
    uri_column: str = DEFAULT_URI_COLUMN,
) -> dict[str, Any]:
    """Genome-wide helper: parallel chroms + flare_vcfs for TractorMixGenome."""
    if not scan_chroms:
        raise SystemExit("scan_chroms must be non-empty")
    missing = [c for c in [*scan_chroms, *grm_chroms] if c not in table]
    if missing:
        raise SystemExit(f"Chromosomes not in table: {missing}; have {sorted(table)}")

    return {
        "chroms": list(scan_chroms),
        "flare_vcfs": [_uri_for(table, c, uri_column) for c in scan_chroms],
        "grm_vcfs": [_uri_for(table, c, uri_column) for c in grm_chroms],
        "grm_chroms": list(grm_chroms),
        "uri_column": uri_column,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--chrom-tsv", type=Path, help="Local TSV export of aou_lr_chrom")
    src.add_argument(
        "--from-firecloud",
        action="store_true",
        help="Fetch aou_lr_chrom via firecloud API",
    )
    src.add_argument(
        "--entities-json",
        type=Path,
        help="Path to saved firecloud get_entities JSON",
    )
    p.add_argument("--namespace", default=DEFAULT_NAMESPACE)
    p.add_argument("--workspace", default=DEFAULT_WORKSPACE)
    p.add_argument("--entity-type", default=DEFAULT_ENTITY_TYPE)
    p.add_argument(
        "--scan-chrom",
        default=None,
        help="Single scan chrom for TractorMixPilot (default chr22 if neither "
        "--scan-chroms nor --autosomes)",
    )
    p.add_argument(
        "--scan-chroms",
        nargs="+",
        default=None,
        help="Scan chrom list for TractorMixGenome (parallel to flare_vcfs)",
    )
    p.add_argument(
        "--autosomes",
        action="store_true",
        help="Shortcut: --scan-chroms chr1 .. chr22 for TractorMixGenome",
    )
    p.add_argument("--grm-chroms", nargs="+", default=["chr1", "chr22"])
    p.add_argument("--uri-column", default=DEFAULT_URI_COLUMN)
    p.add_argument("--out-json", type=Path)
    args = p.parse_args()

    if args.chrom_tsv:
        table = load_chrom_table_tsv(args.chrom_tsv)
    elif args.entities_json:
        entities = json.loads(args.entities_json.read_text())
        table = load_chrom_table_entities(entities)
    else:
        table = fetch_chrom_table_firecloud(args.namespace, args.workspace, args.entity_type)

    if args.autosomes and args.scan_chroms:
        raise SystemExit("Use only one of --autosomes or --scan-chroms")
    if args.autosomes:
        scan_chroms = list(DEFAULT_AUTOSOMES)
    else:
        scan_chroms = args.scan_chroms

    if scan_chroms is not None:
        payload = resolve_genome_uris(
            table,
            scan_chroms=scan_chroms,
            grm_chroms=args.grm_chroms,
            uri_column=args.uri_column,
        )
    else:
        payload = resolve_uris(
            table,
            scan_chrom=args.scan_chrom or "chr22",
            grm_chroms=args.grm_chroms,
            uri_column=args.uri_column,
        )
    print(json.dumps(payload, indent=2))
    if args.out_json:
        args.out_json.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
