#!/usr/bin/env bash
# Copy bam_to_contig FASTAs from a Terra sample table to Julie's DRC prefix.
#
# Reads final_out_{cel_celp,cyp2d6_cyp2d7,triplet_cds,triplet_utr} and writes:
#
#   gs://prod-drc-broad/longreads/wertz/bam_to_contig/cel_celp/{entity}.fa
#   gs://prod-drc-broad/longreads/wertz/bam_to_contig/cyp2d6_7/{entity}.fa
#   gs://prod-drc-broad/longreads/wertz/bam_to_contig/jiadong_cds/{entity}.fa
#   gs://prod-drc-broad/longreads/wertz/bam_to_contig/jiadong_utr/{entity}.fa
#   gs://prod-drc-broad/longreads/wertz/bam_to_contig/copy_manifest.tsv
#
# Source objects are all named all_contigs_concat.fa, so each copy is renamed
# to the Terra entity id ({research_id}_{cohort}). Empty final_out cells are
# skipped (logged). Run on a Terra Workbench VM — laptop gsutil cannot read
# the workspace bucket.
#
# Usage:
#   ./scripts/bam_to_contig/copy_results_to_wertz.sh
#   ./scripts/bam_to_contig/copy_results_to_wertz.sh --dry-run
#   ./scripts/bam_to_contig/copy_results_to_wertz.sh --table sample-hifi-hg38-all-cohorts.tsv
#   JOBS=16 ./scripts/bam_to_contig/copy_results_to_wertz.sh --job cel_celp --limit 20
#
# Requires: bash, python3, gsutil.

set -euo pipefail

DEST="${DEST:-gs://prod-drc-broad/longreads/wertz/bam_to_contig}"
TABLE="${TABLE:-}"
JOBS="${JOBS:-32}"
DRY_RUN=0
LIMIT=""
ONLY_JOB=""
SKIP_EXISTING=0

usage() {
  sed -n '2,24p' "$0" | sed 's/^# \?//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1 ;;
    --table) TABLE="$2"; shift ;;
    --dest) DEST="$2"; shift ;;
    --job) ONLY_JOB="$2"; shift ;;
    --limit) LIMIT="$2"; shift ;;
    --skip-existing) SKIP_EXISTING=1 ;;
    --jobs) JOBS="$2"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

command -v gsutil >/dev/null || { echo "gsutil not found on PATH." >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found on PATH." >&2; exit 1; }

if [[ -z "$TABLE" ]]; then
  for cand in \
    sample-hifi-hg38-all-cohorts.tsv \
    bam_to_contig/configs/sample-hifi-hg38-all-cohorts.tsv \
    "$(dirname "$0")/../../sample-hifi-hg38-all-cohorts.tsv"
  do
    if [[ -f "$cand" ]]; then
      TABLE="$cand"
      break
    fi
  done
fi
if [[ -z "$TABLE" || ! -f "$TABLE" ]]; then
  echo "Set --table / TABLE to the Terra sample-hifi-hg38-all-cohorts TSV." >&2
  exit 1
fi

DEST="${DEST%/}"
GSUTIL=(gsutil)
if [[ -n "${GOOGLE_PROJECT:-}" ]]; then
  GSUTIL=(gsutil -u "$GOOGLE_PROJECT")
fi

export TABLE DEST JOBS DRY_RUN LIMIT ONLY_JOB SKIP_EXISTING
export GSUTIL_BIN="${GSUTIL[0]}"
export GSUTIL_PROJECT="${GOOGLE_PROJECT:-}"

python3 - <<'PY'
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

TABLE = Path(os.environ["TABLE"])
DEST = os.environ["DEST"].rstrip("/")
JOBS = int(os.environ["JOBS"])
DRY_RUN = os.environ["DRY_RUN"] == "1"
LIMIT = os.environ.get("LIMIT") or ""
ONLY_JOB = os.environ.get("ONLY_JOB") or ""
SKIP_EXISTING = os.environ["SKIP_EXISTING"] == "1"
GSUTIL_BIN = os.environ["GSUTIL_BIN"]
GSUTIL_PROJECT = os.environ.get("GSUTIL_PROJECT") or ""

JOB_MAP = {
    "final_out_cel_celp": "cel_celp",
    "final_out_cyp2d6_cyp2d7": "cyp2d6_7",
    "final_out_triplet_cds": "jiadong_cds",
    "final_out_triplet_utr": "jiadong_utr",
}
ID_COL = "entity:sample-hifi-hg38-all-cohorts_id"

gsutil = [GSUTIL_BIN]
if GSUTIL_PROJECT:
    gsutil += ["-u", GSUTIL_PROJECT]


def gsutil_cp(src: str, dest: str) -> subprocess.CompletedProcess:
    cmd = gsutil + ["-q", "cp"]
    if SKIP_EXISTING:
        cmd.append("-n")
    cmd += [src, dest]
    return subprocess.run(cmd, check=False, capture_output=True, text=True)


with TABLE.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

missing_cols = [c for c in JOB_MAP if c not in (rows[0].keys() if rows else [])]
if missing_cols:
    sys.exit(f"Table is missing columns: {', '.join(missing_cols)}")

pairs = []
skipped = []
for r in rows:
    entity = (r.get(ID_COL) or "").strip()
    if not entity:
        continue
    for col, job in JOB_MAP.items():
        if ONLY_JOB and job != ONLY_JOB and col != ONLY_JOB:
            continue
        src = (r.get(col) or "").strip()
        dest = f"{DEST}/{job}/{entity}.fa"
        if not src:
            skipped.append((entity, job, "empty_final_out"))
            continue
        if not src.startswith("gs://"):
            skipped.append((entity, job, f"not_gcs:{src[:80]}"))
            continue
        pairs.append((src, dest, entity, job))

if LIMIT:
    pairs = pairs[: int(LIMIT)]

print(f"table\t{TABLE}")
print(f"dest\t{DEST}")
print(f"copies\t{len(pairs)}")
print(f"skipped_empty\t{len(skipped)}")
print(f"jobs\t{JOBS}")
print(f"dry_run\t{int(DRY_RUN)}")
if skipped:
    print("missing_entities:")
    for entity, job, why in skipped:
        print(f"  {entity}\t{job}\t{why}")

if DRY_RUN:
    for src, dest, entity, job in pairs[:20]:
        print(f"would_copy\t{entity}\t{job}\t{src}\t{dest}")
    if len(pairs) > 20:
        print(f"... {len(pairs) - 20} more")
    sys.exit(0)

if not pairs:
    sys.exit("No final_out URIs to copy.")

ok = 0
fail = []
with ThreadPoolExecutor(max_workers=JOBS) as pool:
    futs = {pool.submit(gsutil_cp, src, dest): (src, dest, entity, job)
            for src, dest, entity, job in pairs}
    done = 0
    for fut in as_completed(futs):
        src, dest, entity, job = futs[fut]
        done += 1
        proc = fut.result()
        if proc.returncode == 0:
            ok += 1
        else:
            err = (proc.stderr or proc.stdout or "").strip().replace("\n", " ")
            fail.append((entity, job, src, dest, err[:300]))
        if done % 500 == 0 or done == len(pairs):
            print(f"progress\t{done}/{len(pairs)}\tok={ok}\tfail={len(fail)}", flush=True)

print(f"copied\t{ok}")
print(f"failed\t{len(fail)}")
for entity, job, src, dest, err in fail:
    print(f"FAIL\t{entity}\t{job}\t{src}\t{dest}\t{err}")

manifest_body = "entity\tjob\tsrc\tdest\n"
for src, dest, entity, job in pairs:
    manifest_body += f"{entity}\t{job}\t{src}\t{dest}\n"
man_local = Path("/tmp/bam_to_contig_wertz_manifest.tsv")
man_local.write_text(manifest_body)
man_dest = f"{DEST}/copy_manifest.tsv"
man = subprocess.run(gsutil + ["-q", "cp", str(man_local), man_dest], check=False)
if man.returncode == 0:
    print(f"manifest\t{man_dest}")
else:
    print("manifest_upload_failed", man.stderr, file=sys.stderr)

sys.exit(1 if fail else 0)
PY
