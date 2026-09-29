#!/usr/bin/env python3
"""Build dsub --tasks TSV files from a per-sample CSV.

The WGS CRAM stays an --env gs:// URI (never --input) so Cloud Batch does not
localize the Nearline object. The CRAI is --env gs:// too: it lives in the
same requester-pays bucket, and dsub/Cromwell localization of that index is
what was failing MakeMinicram.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _gs(prefix: str, sample_id: str, name: str) -> str:
    return f"{prefix.rstrip('/')}/{sample_id}/{name}"


def eh_rows(csv_path: Path, out_prefix: str) -> tuple[list[str], list[list[str]], list[str], list[list[str]]]:
    mini_header = [
        "--env SAMPLE_ID",
        "--env CRAM",
        "--env CRAI",
        "--input REF_FA",
        "--input REF_FAI",
        "--input CATALOG",
        "--env GCLOUD_PROJECT",
        "--env WORKER_LOG_GCS",
        "--output MINICRAM",
        "--output MINICRAI",
        "--output TRANSFER_STATS",
    ]
    gt_header = [
        "--env SAMPLE_ID",
        "--env SEX",
        "--input MINICRAM",
        "--input MINICRAI",
        "--input REF_FA",
        "--input REF_FAI",
        "--input CATALOG",
        "--output EH_JSON",
        "--output EH_VCF",
    ]
    mini_rows: list[list[str]] = []
    gt_rows: list[list[str]] = []
    with csv_path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            sid = rec["sample_id"]
            proj = rec.get("gcloud_project") or rec.get("GCLOUD_PROJECT") or ""
            mini_rows.append(
                [
                    sid,
                    rec["cram"],
                    rec["crai"],
                    rec["ref_fa"],
                    rec["ref_fai"],
                    rec["catalog"],
                    proj,
                    _gs(out_prefix, sid, f"{sid}.minicram.worker.log"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    _gs(out_prefix, sid, f"{sid}.data_transfer_stats.tsv"),
                ]
            )
            gt_rows.append(
                [
                    sid,
                    rec.get("sex") or "female",
                    _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    rec["ref_fa"],
                    rec["ref_fai"],
                    rec["catalog"],
                    _gs(out_prefix, sid, f"{sid}.EH.json"),
                    _gs(out_prefix, sid, f"{sid}.EH.vcf"),
                ]
            )
    return mini_header, mini_rows, gt_header, gt_rows


def locityper_rows(csv_path: Path, out_prefix: str) -> tuple[list[str], list[list[str]], list[str], list[list[str]]]:
    mini_header = [
        "--env SAMPLE_ID",
        "--env CRAM",
        "--env CRAI",
        "--input REF_FA",
        "--input BED",
        "--env GCLOUD_PROJECT",
        "--env WORKER_LOG_GCS",
        "--output MINICRAM",
        "--output MINICRAI",
        "--output TRANSFER_STATS",
    ]
    gt_header = [
        "--env SAMPLE_ID",
        "--input MINICRAM",
        "--input MINICRAI",
        "--input REF_FA",
        "--input REF_FAI",
        "--input COUNTS_JF",
        "--input BED",
        "--input DB_TAR",
        "--output SUMMARY_CSV",
        "--output RESULTS_TAR",
    ]
    mini_rows: list[list[str]] = []
    gt_rows: list[list[str]] = []
    with csv_path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            sid = rec["sample_id"]
            proj = rec.get("gcloud_project") or rec.get("GCLOUD_PROJECT") or ""
            ref_fa = rec.get("ref_fa") or rec["ref_fa_uncompressed"]
            ref_fai = rec.get("ref_fai") or rec["ref_fai_uncompressed"]
            db = rec.get("db_tar") or rec["locityper_db_tar_gz"]
            mini_rows.append(
                [
                    sid,
                    rec["cram"],
                    rec["crai"],
                    ref_fa,
                    rec["bed"],
                    proj,
                    _gs(out_prefix, sid, f"{sid}.minicram.worker.log"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    _gs(out_prefix, sid, f"{sid}.data_transfer_stats.tsv"),
                ]
            )
            gt_rows.append(
                [
                    sid,
                    _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    ref_fa,
                    ref_fai,
                    rec["counts_jf"],
                    rec["bed"],
                    db,
                    _gs(out_prefix, sid, f"{sid}.gts.filtered.csv"),
                    _gs(out_prefix, sid, f"{sid}.locityper.tar.gz"),
                ]
            )
    return mini_header, mini_rows, gt_header, gt_rows


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pipeline", choices=("expansion_hunter", "locityper"), required=True)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--out-prefix", required=True, help="gs://bucket/prefix for per-sample outputs")
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    if args.pipeline == "expansion_hunter":
        mh, mr, gh, gr = eh_rows(args.csv, args.out_prefix)
    else:
        mh, mr, gh, gr = locityper_rows(args.csv, args.out_prefix)
    write_tsv(args.outdir / "minicram.tasks.tsv", mh, mr)
    write_tsv(args.outdir / "genotype.tasks.tsv", gh, gr)
    print(f"wrote {len(mr)} minicram tasks and {len(gr)} genotype tasks under {args.outdir}")


if __name__ == "__main__":
    main()
