#!/usr/bin/env bash
# Collect every monitoring.log under a Cromwell/Terra workflow root into a tarball.
#
# Terra tasks in this workspace localize cromwell_monitoring_script.sh, so each
# task execution directory holds a monitoring.log with the #CPU / Total Memory /
# Total Disk header and a 10-second CPU/memory/disk/IO trace. The log carries no
# task name, so the GCS path is the only key: the tarball mirrors the path under
# the workflow root (call-<Task>/shard-N/attempt-K/monitoring.log, including any
# subworkflow levels) and ships a manifest.tsv keyed the same way.
#
# Writes <OUT>.tar.gz containing:
#   <OUT>/manifest.tsv                    relpath  call  shard  attempt  subworkflow  size_bytes  gs_uri
#   <OUT>/object_index.tsv                relpath  size_bytes  created_utc  -- every object
#                                         under the run root. monitoring.log covers only the
#                                         task command, so creation times of gcs_localization.sh
#                                         (job start), rc (command end) and the delocalized
#                                         outputs are what recover localization/delocalization
#                                         time and output sizes. Skip with --no-index.
#   <OUT>/run_root.txt                    the gs:// root these logs came from
#   <OUT>/call-<Task>/shard-N/[attempt-K/]monitoring.log
#
# Usage:
#   ./scripts/collect_terra_monitoring_logs.sh gs://fc-secure-.../FelixGenome/<workflow-id>
#   FELIX_RUN_ROOT=gs://... ./scripts/collect_terra_monitoring_logs.sh
#   WORKSPACE_BUCKET=gs://fc-secure-... FELIX_SUBMISSION_ID=<sid> FELIX_WORKFLOW_ID=<wid> \
#     WORKFLOW_NAME=FelixGenome ./scripts/collect_terra_monitoring_logs.sh
#   OUT=felix_monitoring JOBS=16 ./scripts/collect_terra_monitoring_logs.sh gs://...
#   ./scripts/collect_terra_monitoring_logs.sh --dry-run gs://...    # list only, no copy
#   ./scripts/collect_terra_monitoring_logs.sh --no-index gs://...   # skip the object index
#
# Set GOOGLE_PROJECT to bill a requester-pays bucket (adds gsutil -u).
#
# Requires: bash, gsutil, awk, xargs, tar.

set -euo pipefail

RUN_ROOT="${FELIX_RUN_ROOT:-}"
OUT="${OUT:-terra_monitoring}"
JOBS="${JOBS:-16}"
WORKFLOW_NAME="${WORKFLOW_NAME:-FelixGenome}"
DRY_RUN=0
INDEX="${INDEX:-1}"

usage() { sed -n '2,33p' "$0" | sed 's/^# \?//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1 ;;
    --no-index) INDEX=0 ;;
    -j|--jobs) JOBS="$2"; shift ;;
    -o|--out) OUT="$2"; shift ;;
    -h|--help) usage; exit 0 ;;
    gs://*) RUN_ROOT="$1" ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "${RUN_ROOT}" ]]; then
  if [[ -n "${WORKSPACE_BUCKET:-}" && -n "${FELIX_SUBMISSION_ID:-}" && -n "${FELIX_WORKFLOW_ID:-}" ]]; then
    RUN_ROOT="${WORKSPACE_BUCKET%/}/submissions/${FELIX_SUBMISSION_ID}/${WORKFLOW_NAME}/${FELIX_WORKFLOW_ID}"
  else
    echo "No workflow root. Pass gs://... or set FELIX_RUN_ROOT (or WORKSPACE_BUCKET + FELIX_SUBMISSION_ID + FELIX_WORKFLOW_ID)." >&2
    exit 1
  fi
fi
RUN_ROOT="${RUN_ROOT%/}"

# Plain function, not an array: bash 3.2 (macOS) errors on an empty "${arr[@]}"
# under `set -u`.
gs() {
  if [[ -n "${GOOGLE_PROJECT:-}" ]]; then gsutil -u "${GOOGLE_PROJECT}" "$@"; else gsutil "$@"; fi
}

# Explicit template: BSD mktemp (macOS) ignores TMPDIR for a bare `mktemp -d`.
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/terra_monitoring.XXXXXX")"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Listing monitoring.log under ${RUN_ROOT} ..." >&2
# -l gives "size  date  uri" so the manifest records sizes without a second pass.
gs ls -l "${RUN_ROOT}/**monitoring.log" 2>/dev/null \
  | awk -v root="${RUN_ROOT}/" '
      $NF ~ /^gs:\/\// && $NF ~ /monitoring\.log$/ {
        uri = $NF; size = $1
        rel = uri; sub("^" root, "", rel)
        print rel "\t" size "\t" uri
      }' \
  | sort > "${WORKDIR}/index.tsv" || true

N="$(wc -l < "${WORKDIR}/index.tsv" | tr -d ' ')"
if [[ "${N}" -eq 0 ]]; then
  echo "No monitoring.log objects found under ${RUN_ROOT}" >&2
  exit 1
fi
echo "Found ${N} monitoring.log files." >&2

# relpath -> call / shard / attempt / subworkflow. The last call-* segment is the
# task that actually ran; anything before it is subworkflow nesting.
{
  printf 'relpath\tcall\tshard\tattempt\tsubworkflow\tsize_bytes\tgs_uri\n'
  awk -F'\t' '{
    n = split($1, p, "/")
    call = ""; shard = "NA"; attempt = "1"; sub_wf = "NA"; last_call = 0
    for (i = 1; i <= n; i++) {
      if (p[i] ~ /^call-/) { call = substr(p[i], 6); last_call = i }
      else if (p[i] ~ /^shard-/ && i > last_call) { shard = substr(p[i], 7) }
      else if (p[i] ~ /^attempt-/) { attempt = substr(p[i], 9) }
    }
    for (i = 1; i < last_call; i++) if (p[i] ~ /^call-/) sub_wf = (sub_wf == "NA" ? substr(p[i], 6) : sub_wf "/" substr(p[i], 6))
    if (call == "") call = "UNKNOWN"
    print $1 "\t" call "\t" shard "\t" attempt "\t" sub_wf "\t" $2 "\t" $3
  }' "${WORKDIR}/index.tsv"
} > "${WORKDIR}/manifest.tsv"

echo "Tasks:" >&2
awk -F'\t' 'NR>1 {c[$2]++} END {for (k in c) printf "  %-28s %5d\n", k, c[k]}' "${WORKDIR}/manifest.tsv" | sort >&2

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "--dry-run: manifest below, nothing copied." >&2
  cat "${WORKDIR}/manifest.tsv"
  exit 0
fi

STAGE="${WORKDIR}/${OUT}"
mkdir -p "${STAGE}"
cp "${WORKDIR}/manifest.tsv" "${STAGE}/manifest.tsv"
printf '%s\n' "${RUN_ROOT}" > "${STAGE}/run_root.txt"

# One recursive listing of the whole run root. `gsutil ls -l` prints
# "<size> <ISO-8601 create time> <uri>" per object (and a trailing TOTAL line,
# which has no gs:// in the last field and is dropped). Objects created by the
# worker are timestamped at upload, so within a task directory the spread of
# creation times brackets the command and its delocalization.
if [[ "${INDEX}" -eq 1 ]]; then
  echo "Indexing all objects under ${RUN_ROOT} ..." >&2
  {
    printf 'relpath\tsize_bytes\tcreated_utc\n'
    gs ls -l "${RUN_ROOT}/**" 2>/dev/null \
      | awk -v root="${RUN_ROOT}/" '
          NF >= 3 && $NF ~ /^gs:\/\// && $1 ~ /^[0-9]+$/ {
            uri = $NF; rel = uri; sub("^" root, "", rel)
            print rel "\t" $1 "\t" $2
          }' \
      | sort
  } > "${STAGE}/object_index.tsv"
  echo "Indexed $(( $(wc -l < "${STAGE}/object_index.tsv") - 1 )) objects." >&2
fi

# One gsutil per file, JOBS at a time. gsutil -m cp -I would be one process but
# flattens the names, and the path is the only thing identifying the task.
export STAGE
[[ -n "${GOOGLE_PROJECT:-}" ]] && export GOOGLE_PROJECT
fetch_one() {
  rel="$1"; uri="$2"
  dest="${STAGE}/${rel}"
  mkdir -p "$(dirname "${dest}")"
  if [[ -n "${GOOGLE_PROJECT:-}" ]]; then
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
  else
    gsutil -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
  fi
}
export -f fetch_one

# -n2 (not -I) so relpath and URI arrive as two separate arguments; neither
# Cromwell paths nor gs:// URIs contain whitespace.
cut -f1,3 "${WORKDIR}/index.tsv" | tr '\t' '\n' \
  | xargs -P "${JOBS}" -n2 bash -c 'fetch_one "$1" "$2"' _

GOT="$(find "${STAGE}" -name monitoring.log | wc -l | tr -d ' ')"
echo "Downloaded ${GOT}/${N} logs." >&2

tar -czf "${OUT}.tar.gz" -C "${WORKDIR}" "${OUT}"
echo "Wrote ${OUT}.tar.gz ($(du -h "${OUT}.tar.gz" | cut -f1), ${GOT} logs)" >&2
