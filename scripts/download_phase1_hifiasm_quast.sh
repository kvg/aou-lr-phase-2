#!/usr/bin/env bash
# Download Phase 1 hifiasm QUAST reports and melt them into one TSV.
#
# Source objects are all named report_map.txt; the research_id is the parent
# directory. This writes one file per sample plus a long table:
#
#   <research_id>.report_map.txt
#   hifiasm_quast.tsv    research_id  assembly  metric  value
#
# Usage:
#   export GOOGLE_PROJECT=your-requester-pays-billing-project
#   ./scripts/download_phase1_hifiasm_quast.sh
#   ./scripts/download_phase1_hifiasm_quast.sh --dry-run
#   JOBS=32 OUT=/tmp/quast ./scripts/download_phase1_hifiasm_quast.sh
#
# Requires: bash, gsutil, awk, xargs.

set -euo pipefail

SRC_PREFIX="${SRC_PREFIX:-gs://vwb-aou-datasets-controlled/pooled/longreads/v7_base/hifiasm}"
OUT="${OUT:-hifiasm_quast}"
JOBS="${JOBS:-16}"
DRY_RUN=0

usage() {
  sed -n '2,18p' "$0" | sed 's/^# \?//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1 ;;
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

GSUTIL=(gsutil -u "$GOOGLE_PROJECT")
LIST_GLOB="${SRC_PREFIX%/}/*/report_map.txt"

URI_LIST="$(mktemp)"
FAIL_LIST="$(mktemp)"
trap 'rm -f "$URI_LIST" "$FAIL_LIST"' EXIT

echo "Listing ${LIST_GLOB}"
"${GSUTIL[@]}" ls "$LIST_GLOB" > "$URI_LIST"
n_uri="$(wc -l < "$URI_LIST" | tr -d " ")"
if [[ "$n_uri" -eq 0 ]]; then
  echo "No report_map.txt objects found." >&2
  exit 1
fi
echo "Found ${n_uri} reports"

if [[ "$DRY_RUN" -eq 1 ]]; then
  sed 's#.*/hifiasm/\([^/]*\)/report_map.txt#\1#' "$URI_LIST" | sort -u | head
  echo "... ${n_uri} sample directories"
  exit 0
fi

mkdir -p "$OUT"

copy_one() {
  local uri="$1"
  local sid dest
  sid=$(basename "$(dirname "$uri")")
  dest="${OUT}/${sid}.report_map.txt"
  if [[ -s "$dest" ]]; then
    return 0
  fi
  gsutil -q -u "$GOOGLE_PROJECT" cp "$uri" "$dest" || echo "$uri" >> "$FAIL_LIST"
}
export -f copy_one
export OUT FAIL_LIST GOOGLE_PROJECT

xargs -P "$JOBS" -n 1 bash -c 'copy_one "$1"' _ < "$URI_LIST"

n_fail="$(wc -l < "$FAIL_LIST" | tr -d " ")"
if [[ "$n_fail" -gt 0 ]]; then
  echo "Failed to copy ${n_fail} objects:" >&2
  cat "$FAIL_LIST" >&2
  exit 1
fi

n_files="$(find "$OUT" -maxdepth 1 -name "*.report_map.txt" -type f | wc -l | tr -d " ")"
echo "Wrote ${n_files} files under ${OUT}"

COMBINED="${OUT}/hifiasm_quast.tsv"
{
  printf "research_id\tassembly\tmetric\tvalue\n"
  for report in "$OUT"/*.report_map.txt; do
    sid="$(basename "$report" .report_map.txt)"
    awk -F '\t' -v sid="$sid" '
      NR == 1 {
        for (i = 2; i <= NF; i++) {
          name = $i
          sub("^" sid "\\.", "", name)
          sub("^bp\\.", "", name)
          sub("\\.p_ctg$", "", name)
          if (name == "" || name == "p_ctg") name = "p_ctg"
          col[i] = name
        }
        next
      }
      NF < 2 { next }
      {
        metric = $1
        for (i = 2; i <= NF; i++) {
          printf "%s\t%s\t%s\t%s\n", sid, col[i], metric, $i
        }
      }
    ' "$report"
  done
} > "$COMBINED"

echo "Wrote ${COMBINED}"
echo "Rows: $(($(wc -l < "$COMBINED") - 1))"
echo "Unique research_ids: $(awk -F '\t' 'NR>1 {print $1}' "$COMBINED" | sort -u | wc -l | tr -d " ")"
