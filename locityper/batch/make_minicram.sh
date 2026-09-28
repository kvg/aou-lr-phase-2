#!/usr/bin/env bash
# Cloud Batch / dsub: Nearline-efficient minicram for Locityper (print_reads).
# CRAM must be a gs:// URI (--env), not a localized --input.
set -euo pipefail
date

: "${SAMPLE_ID:?}"
: "${CRAM:?}"
: "${CRAI:?}"
: "${REF_FA:?}"
: "${BED:?}"
: "${MINICRAM:?}"
: "${MINICRAI:?}"
: "${TRANSFER_STATS:?}"

WINDOW_GRAB="${WINDOW_GRAB:-3000}"
# GRCh38 Locityper background interval (Isaac). Empty skips the extra grab.
BG_REGION_BED="${BG_REGION_BED:-chr17	72062001	76562000}"
MAX_RETRY="${MAX_RETRY:-3}"
WAIT_TIME="${WAIT_TIME:-30}"

awk -v w="${WINDOW_GRAB}" 'BEGIN { OFS="\t" } {
  s = $2 - w
  if (s < 0) s = 0
  print $1, s, $3 + w
}' "${BED}" > intervals.bed

if [[ -n "${BG_REGION_BED}" ]]; then
  echo -e "${BG_REGION_BED}" | awk -v w="${WINDOW_GRAB}" 'BEGIN { OFS="\t" } {
    s = $2 - w
    if (s < 0) s = 0
    print $1, s, $3 + w
  }' >> intervals.bed
fi

echo "minicram intervals:"
cat intervals.bed
wc -l intervals.bed

PROJ="${GCLOUD_PROJECT:-${GOOGLE_CLOUD_PROJECT:-}}"
if [[ -n "${PROJ}" ]]; then
  echo "requester-pays project: ${PROJ}"
fi

tmp_cram="${TMPDIR:-/tmp}/${SAMPLE_ID}.minicram.cram"

run_print_reads() {
  local -a cmd=(
    python3 -m str_analysis.print_reads
    -R "${REF_FA}"
    --read-index "${CRAI}"
    -L intervals.bed
    --padding 0
    -o "${tmp_cram}"
    --verbose
    --output-data-transfer-stats
  )
  if [[ -n "${PROJ}" ]]; then
    cmd+=(--gcloud-project "${PROJ}")
  fi
  cmd+=("${CRAM}")
  "${cmd[@]}"
}

CURR=0
until (( CURR == MAX_RETRY )); do
  if run_print_reads; then
    break
  fi
  echo "print_reads failed; retry $((CURR + 1))/${MAX_RETRY}" >&2
  rm -f "${tmp_cram}" "${tmp_cram}.crai" *.data_transfer_stats.tsv
  CURR=$((CURR + 1))
  if (( CURR == MAX_RETRY )); then
    echo "print_reads failed after ${MAX_RETRY} attempts" >&2
    exit 1
  fi
  sleep "${WAIT_TIME}"
done

if [[ ! -s "${tmp_cram}" || ! -s "${tmp_cram}.crai" ]]; then
  echo "print_reads wrote no CRAM/CRAI" >&2
  exit 1
fi

mkdir -p "$(dirname "${MINICRAM}")" "$(dirname "${MINICRAI}")" "$(dirname "${TRANSFER_STATS}")"
cp -f "${tmp_cram}" "${MINICRAM}"
cp -f "${tmp_cram}.crai" "${MINICRAI}"
stats=$(ls -1 *.data_transfer_stats.tsv | head -n 1)
cp -f "${stats}" "${TRANSFER_STATS}"
ls -lh "${MINICRAM}" "${MINICRAI}" "${TRANSFER_STATS}"
date
