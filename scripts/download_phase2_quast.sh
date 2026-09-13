#!/usr/bin/env bash
# Download every Phase 2 hifiasm QUAST summary from the v9 long-read manifest.
#
# URIs come from column assembly_quast_report_summary. People can have more
# than one report (v8 vs v9, PacBio vs ONT, multiple GCs). This keeps them
# all; join later against covariates. Writes:
#
#   <release>/<GC>/<platform>/assembly/<id>.quast-summary.txt
#   phase2_quast.tsv           research_id  release  gc  platform  uri  assembly  metric  value
#   phase2_quast.sources.tsv   research_id  release  gc  platform  relpath  uri
#
# Copies with gsutil -m per unique .../assembly/ prefix (one process, many
# objects). JOBS is how many prefixes to copy at once.
#
# Usage:
#   export GOOGLE_PROJECT=your-requester-pays-billing-project
#   ./scripts/download_phase2_quast.sh
#   ./scripts/download_phase2_quast.sh --dry-run
#   MANIFEST=/path/to/manifest.tsv JOBS=8 ./scripts/download_phase2_quast.sh
#
# Requires: bash, gsutil, awk, xargs, sort, python3.

set -euo pipefail

MANIFEST="${MANIFEST:-}"
QUAST_COL="${QUAST_COL:-assembly_quast_report_summary}"
OUT="${OUT:-phase2_quast}"
JOBS="${JOBS:-8}"
DRY_RUN=0

usage() {
  sed -n '2,22p' "$0" | sed 's/^# \?//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1 ;;
    --manifest) MANIFEST="$2"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "${GOOGLE_PROJECT:-}" ]]; then
  echo "Set GOOGLE_PROJECT to a billing project that can read the requester-pays bucket." >&2
  exit 1
fi
command -v gsutil >/dev/null || { echo "gsutil not found on PATH." >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found on PATH." >&2; exit 1; }

if [[ -z "$MANIFEST" ]]; then
  for cand in \
    workspace/vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv \
    gs://vwb-aou-datasets-controlled-v9/v9/wgs/long_read/manifest.tsv
  do
    if [[ "$cand" == gs://* ]] || [[ -f "$cand" ]]; then
      MANIFEST="$cand"
      break
    fi
  done
fi
if [[ -z "$MANIFEST" ]]; then
  echo "Set MANIFEST to the v9 long-read manifest TSV (local path or gs:// URI)." >&2
  exit 1
fi

GSUTIL=(gsutil -u "$GOOGLE_PROJECT")
SRC_LIST="$(mktemp)"
PREFIX_LIST="$(mktemp)"
FAIL_LIST="$(mktemp)"
trap 'rm -f "$SRC_LIST" "$PREFIX_LIST" "$FAIL_LIST"' EXIT

read_manifest() {
  if [[ "$MANIFEST" == gs://* ]]; then
    "${GSUTIL[@]}" cat "$MANIFEST"
  else
    cat "$MANIFEST"
  fi
}

echo "Reading ${QUAST_COL} from ${MANIFEST}"
read_manifest | awk -F '\t' -v colname="$QUAST_COL" -v OFS='\t' '
  NR == 1 {
    sub(/\r$/, "", $0)
    for (i = 1; i <= NF; i++) {
      name = $i
      sub(/\r$/, "", name)
      if (name == colname) col = i
    }
    if (!col) {
      printf "Missing column %s\n", colname > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    sub(/\r$/, "", $col)
    uri = $col
    if (uri == "" || uri == "NA" || uri == "na") next
    n = split(uri, a, "/")
    sid = a[n]
    sub(/\.quast-summary\.txt$/, "", sid)
    if (sid == "" || sid == a[n]) next
    start = 0
    for (i = 1; i <= n; i++) if (a[i] == "longreads") start = i + 1
    if (start == 0) start = n - 4
    if (start < 1) start = 1
    release = (start <= n) ? a[start] : ""
    gc = (start + 1 <= n) ? a[start + 1] : ""
    platform = (start + 2 <= n) ? a[start + 2] : ""
    rel = release "/" gc "/" platform "/assembly/" sid ".quast-summary.txt"
    print sid, release, gc, platform, rel, uri
  }
' | sort -u > "$SRC_LIST"

if [[ ! -s "$SRC_LIST" ]]; then
  echo "No QUAST URIs found in ${QUAST_COL}." >&2
  exit 1
fi

awk -F '\t' '{
  uri = $6
  n = split(uri, a, "/")
  prefix = a[1]
  for (i = 2; i <= n - 2; i++) prefix = prefix "/" a[i]
  print prefix
}' "$SRC_LIST" | sort -u > "$PREFIX_LIST"

n_uri="$(wc -l < "$SRC_LIST" | tr -d " ")"
n_people="$(awk -F '\t' '{ print $1 }' "$SRC_LIST" | sort -u | wc -l | tr -d " ")"
n_prefix="$(wc -l < "$PREFIX_LIST" | tr -d " ")"
echo "Found ${n_uri} QUAST URIs for ${n_people} research_ids across ${n_prefix} prefixes (keeping all reports)"

if [[ "$DRY_RUN" -eq 1 ]]; then
  head "$PREFIX_LIST"
  echo "... ${n_prefix} prefixes, ${n_uri} URIs"
  exit 0
fi

mkdir -p "$OUT"
{
  printf "research_id\trelease\tgc\tplatform\trelpath\turi\n"
  cat "$SRC_LIST"
} > "${OUT}/phase2_quast.sources.tsv"

copy_prefix() {
  local prefix="$1"
  local rel dest
  if [[ "$prefix" == *"/longreads/"* ]]; then
    rel="${prefix##*/longreads/}"
  else
    rel=$(basename "$prefix")
  fi
  dest="${OUT}/${rel}"
  mkdir -p "$dest"
  gsutil -m -q -u "$GOOGLE_PROJECT" cp -n "${prefix}/*/*.quast-summary.txt" "${dest}/" \
    || echo "$prefix" >> "$FAIL_LIST"
}
export -f copy_prefix
export OUT FAIL_LIST GOOGLE_PROJECT

echo "Copying with gsutil -m (${JOBS} prefixes at a time)"
xargs -P "$JOBS" -n 1 bash -c 'copy_prefix "$1"' _ < "$PREFIX_LIST"

n_fail="$(wc -l < "$FAIL_LIST" | tr -d " ")"
if [[ "$n_fail" -gt 0 ]]; then
  echo "Failed to copy ${n_fail} prefixes:" >&2
  cat "$FAIL_LIST" >&2
  exit 1
fi

n_files="$(find "$OUT" -name "*.quast-summary.txt" -type f | wc -l | tr -d " ")"
echo "Wrote ${n_files} files under ${OUT}"

COMBINED="${OUT}/phase2_quast.tsv"
python3 - "$OUT" "$COMBINED" << 'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])
combined = Path(sys.argv[2])
sources = out / "phase2_quast.sources.tsv"
missing = []
n_rows = 0
people = set()

with sources.open() as fh, combined.open("w", newline="") as out_fh:
    reader = csv.DictReader(fh, delimiter="\t")
    writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["research_id", "release", "gc", "platform", "uri", "assembly", "metric", "value"]
    )
    for row in reader:
        report = out / row["relpath"]
        if not report.is_file():
            missing.append(row["relpath"])
            continue
        sid = row["research_id"]
        lines = report.read_text().splitlines()
        if not lines:
            missing.append(row["relpath"])
            continue
        header = lines[0].split("\t")
        assemblies = []
        for name in header[1:]:
            name = name[len(sid) + 1 :] if name.startswith(sid + ".") else name
            if name.startswith("bp."):
                name = name[3:]
            if name.endswith(".p_ctg"):
                name = name[: -len(".p_ctg")]
            if name in {"", "p_ctg"}:
                name = "p_ctg"
            assemblies.append(name)
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            metric = parts[0]
            for assembly, value in zip(assemblies, parts[1:]):
                writer.writerow(
                    [
                        sid,
                        row["release"],
                        row["gc"],
                        row["platform"],
                        row["uri"],
                        assembly,
                        metric,
                        value,
                    ]
                )
                n_rows += 1
        people.add(sid)

if missing:
    sys.stderr.write(f"Missing {len(missing)} reports after copy. First: {missing[0]}\n")
    sys.exit(1)
print(f"Wrote {combined}")
print(f"Rows: {n_rows}")
print(f"Unique research_ids: {len(people)}")
PY

echo "Wrote ${OUT}/phase2_quast.sources.tsv"
