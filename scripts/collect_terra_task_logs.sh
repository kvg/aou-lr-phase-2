#!/usr/bin/env bash
# Collect stdout / stderr / rc from every task execution under a Terra workflow
# root, and summarise how far the run has got.
#
# Cromwell writes `rc` only when a task finishes, so the object listing alone
# answers "what is still running": a shard with stdout but no rc is in flight.
# That summary is printed before anything is downloaded, and shipped as
# progress.tsv. The logs themselves mirror the path under the workflow root,
# because the files carry no task name.
#
# Writes <OUT>.tar.gz containing:
#   <OUT>/progress.tsv     call  shards  finished  running  failed  worst_rc
#   <OUT>/manifest.tsv     relpath  kind  call  shard  attempt  subworkflow  size_bytes  gs_uri
#   <OUT>/run_root.txt     the gs:// root these logs came from
#   <OUT>/call-<Task>/shard-N/[attempt-K/]{stdout,stderr,rc}
#
# Usage:
#   ./scripts/collect_terra_task_logs.sh gs://fc-secure-.../FelixGenome/<workflow-id>
#   FELIX_RUN_ROOT=gs://... ./scripts/collect_terra_task_logs.sh
#   WORKSPACE_BUCKET=gs://fc-secure-... FELIX_SUBMISSION_ID=<sid> FELIX_WORKFLOW_ID=<wid> \
#     WORKFLOW_NAME=FelixGenome ./scripts/collect_terra_task_logs.sh
#   ./scripts/collect_terra_task_logs.sh --progress-only gs://...   # listing only, no download
#   ./scripts/collect_terra_task_logs.sh --tail-bytes 20000 gs://.. # last 20 kB of each log
#   ./scripts/collect_terra_task_logs.sh --calls RuLoci,RuStep2 gs://...
#   OUT=felix_logs JOBS=24 ./scripts/collect_terra_task_logs.sh gs://...
#
# Set GOOGLE_PROJECT to bill a requester-pays bucket (adds gsutil -u).
#
# Requires: bash, gsutil, awk, xargs, tar.

set -euo pipefail

RUN_ROOT="${FELIX_RUN_ROOT:-}"
OUT="${OUT:-terra_task_logs}"
JOBS="${JOBS:-16}"
WORKFLOW_NAME="${WORKFLOW_NAME:-FelixGenome}"
TAIL_BYTES=0
CALLS=""
PROGRESS_ONLY=0
WARN_BYTES=$((2 * 1024 * 1024 * 1024))   # warn past 2 GiB of logs

usage() { sed -n '2,33p' "$0" | sed 's/^# \?//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--progress-only) PROGRESS_ONLY=1 ;;
    -t|--tail-bytes) TAIL_BYTES="$2"; shift ;;
    -c|--calls) CALLS="$2"; shift ;;
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
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/terra_task_logs.XXXXXX")"
trap 'rm -rf "${WORKDIR}"' EXIT

echo "Listing stdout / stderr / rc under ${RUN_ROOT} ..." >&2
# One ls for all three names; -l gives "size date uri" so sizes need no second pass.
gs ls -l "${RUN_ROOT}/**/stdout" "${RUN_ROOT}/**/stderr" "${RUN_ROOT}/**/rc" 2>/dev/null \
  | awk -v root="${RUN_ROOT}/" '
      $NF ~ /^gs:\/\// {
        uri = $NF; size = $1
        n = split(uri, q, "/"); kind = q[n]
        if (kind != "stdout" && kind != "stderr" && kind != "rc") next
        rel = uri; sub("^" root, "", rel)
        print rel "\t" kind "\t" size "\t" uri
      }' \
  | sort > "${WORKDIR}/all.tsv" || true

if [[ ! -s "${WORKDIR}/all.tsv" ]]; then
  echo "No stdout/stderr/rc objects found under ${RUN_ROOT}" >&2
  exit 1
fi

# relpath -> call / shard / attempt / subworkflow. The last call-* segment is the
# task that actually ran; anything before it is subworkflow nesting.
awk -F'\t' '{
  n = split($1, p, "/")
  call = ""; shard = "NA"; attempt = "1"; sub_wf = "NA"; last_call = 0
  for (i = 1; i <= n; i++) {
    if (p[i] ~ /^call-/) { call = substr(p[i], 6); last_call = i }
    else if (p[i] ~ /^attempt-/) { attempt = substr(p[i], 9) }
  }
  # A nested scatter gives call-Step2/shard-<chrom>/shard-<pheno>/, so keep every
  # shard index after the last call- segment joined, not just the innermost.
  for (i = last_call + 1; i <= n; i++)
    if (p[i] ~ /^shard-/) shard = (shard == "NA" ? substr(p[i], 7) : shard "." substr(p[i], 7))
  for (i = 1; i < last_call; i++) if (p[i] ~ /^call-/) sub_wf = (sub_wf == "NA" ? substr(p[i], 6) : sub_wf "/" substr(p[i], 6))
  if (call == "") call = "UNKNOWN"
  # The execution directory identifies one task execution: the relpath without
  # the trailing filename and without the attempt segment.
  execdir = ""
  for (i = 1; i <= n - 1; i++) if (p[i] !~ /^attempt-/) execdir = (execdir == "" ? p[i] : execdir "/" p[i])
  print $1 "\t" $2 "\t" call "\t" shard "\t" attempt "\t" sub_wf "\t" $3 "\t" $4 "\t" execdir
}' "${WORKDIR}/all.tsv" > "${WORKDIR}/parsed.tsv"

if [[ -n "${CALLS}" ]]; then
  awk -F'\t' -v want="${CALLS}" '
    BEGIN { n = split(want, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
    keep[$3]' "${WORKDIR}/parsed.tsv" > "${WORKDIR}/parsed.filtered.tsv"
  mv "${WORKDIR}/parsed.filtered.tsv" "${WORKDIR}/parsed.tsv"
  if [[ ! -s "${WORKDIR}/parsed.tsv" ]]; then
    echo "No logs for calls: ${CALLS}" >&2
    exit 1
  fi
fi

{
  printf 'relpath\tkind\tcall\tshard\tattempt\tsubworkflow\tsize_bytes\tgs_uri\texec_dir\n'
  cat "${WORKDIR}/parsed.tsv"
} > "${WORKDIR}/manifest.tsv"

# Progress, one row per SHARD (attempts collapsed). A shard counts as finished
# when its highest attempt has an rc; stdout with no rc on that attempt means it
# is in flight. `expected` is how many shards the task will have in the end,
# derived from the scatter widths the run itself reveals: FitFelixNull runs once
# per phenotype and PackFelixla once per chromosome, so the step-2 scatters are
# their product. That is the difference between "87 done" and "87 of 1,914".
awk -F'\t' '
  { call = $3; attempt = $5 + 0
    key = $9                       # execution directory: unique per task execution
    if (attempt > top[key]) top[key] = attempt
    if ($2 == "rc") rcseen[key "\t" attempt] = 1
    attempts[key "\t" attempt] = 1
    name[key] = call }
  END {
    for (k in name) { c = name[k]; shards[c]++
      if ((k "\t" top[k]) in rcseen) fin[c]++; else run[c]++ }
    # key is now one field, so the attempt is f[2] and the call comes from name[].
    for (a in attempts) { split(a, f, "\t"); if (f[2] + 0 > 1) retried[name[f[1]]]++ }

    P = shards["Null"] + 0; C = shards["Pack"] + 0
    for (c in shards) {
      want = 0
      if (c == "Step2" || c == "RuStep2" || c == "BaStep2") want = C * P
      else if (c == "RuLoci" || c == "RuVcf" || c == "BaVcf" || c == "Pack") want = C
      else if (c == "Null" || c == "Concat" || c == "RuConcat" || c == "BaConcat" || c == "ArmCompare") want = P
      if (want < shards[c]) want = shards[c]
      printf "%s\t%d\t%d\t%d\t%d\t%d\n", c, want, shards[c], fin[c] + 0, run[c] + 0, retried[c] + 0
    }
  }' "${WORKDIR}/parsed.tsv" | sort > "${WORKDIR}/progress.body"

{
  printf 'call\texpected\tstarted\tfinished\trunning\tretried\n'
  cat "${WORKDIR}/progress.body"
} > "${WORKDIR}/progress.tsv"

echo "" >&2
printf '%-14s %9s %9s %9s %9s %8s %7s\n' call expected started finished running retried pct >&2
awk -F'\t' '{printf "%-14s %9d %9d %9d %9d %8d %6.1f%%\n", $1, $2, $3, $4, $5, $6, 100*$4/$2}' \
  "${WORKDIR}/progress.body" >&2
awk -F'\t' '{e+=$2; s+=$3; f+=$4; r+=$5; t+=$6}
  END {printf "%-14s %9d %9d %9d %9d %8d %6.1f%%\n", "TOTAL", e, s, f, r, t, 100*f/e}' \
  "${WORKDIR}/progress.body" >&2

# printf, not print: awk renders a 28 GB total as 2.8e+10 and bash arithmetic
# then refuses it.
TOTAL_BYTES="$(awk -F'\t' '{s+=$7} END {printf "%.0f\n", s+0}' "${WORKDIR}/parsed.tsv")"
TOTAL_MIB="$(awk -F'\t' '{s+=$7} END {printf "%.0f\n", s/1048576}' "${WORKDIR}/parsed.tsv")"
echo "" >&2
echo "$(wc -l < "${WORKDIR}/parsed.tsv" | tr -d ' ') objects, ${TOTAL_MIB} MiB of logs." >&2

if [[ "${PROGRESS_ONLY}" -eq 1 ]]; then
  echo "--progress-only: nothing downloaded." >&2
  cat "${WORKDIR}/progress.tsv"
  exit 0
fi
if [[ "${TAIL_BYTES}" -eq 0 && "${TOTAL_BYTES}" -gt "${WARN_BYTES}" ]]; then
  echo "That is a lot to download. Re-run with --tail-bytes 20000 to keep only the" >&2
  echo "end of each log, or --calls <Task> to narrow it." >&2
fi

STAGE="${WORKDIR}/${OUT}"
mkdir -p "${STAGE}"
cp "${WORKDIR}/manifest.tsv" "${STAGE}/manifest.tsv"
cp "${WORKDIR}/progress.tsv" "${STAGE}/progress.tsv"
printf '%s\n' "${RUN_ROOT}" > "${STAGE}/run_root.txt"

# One gsutil per file, JOBS at a time. gsutil -m cp -I would be one process but
# flattens the names, and the path is the only thing identifying the task.
export STAGE TAIL_BYTES
[[ -n "${GOOGLE_PROJECT:-}" ]] && export GOOGLE_PROJECT
fetch_one() {
  rel="$1"; uri="$2"
  dest="${STAGE}/${rel}"
  mkdir -p "$(dirname "${dest}")"
  # rc is two bytes and is the thing progress is read from, so never truncate it.
  if [[ "${TAIL_BYTES}" -gt 0 && "${rel##*/}" != "rc" ]]; then
    if [[ -n "${GOOGLE_PROJECT:-}" ]]; then
      gsutil -u "${GOOGLE_PROJECT}" cat -r "-${TAIL_BYTES}" "${uri}" > "${dest}" 2>/dev/null \
        || gsutil -u "${GOOGLE_PROJECT}" -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
    else
      gsutil cat -r "-${TAIL_BYTES}" "${uri}" > "${dest}" 2>/dev/null \
        || gsutil -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
    fi
  elif [[ -n "${GOOGLE_PROJECT:-}" ]]; then
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
  else
    gsutil -q cp "${uri}" "${dest}" || echo "FAILED ${uri}" >&2
  fi
}
export -f fetch_one

# -n2 (not -I) so relpath and URI arrive as two separate arguments; neither
# Cromwell paths nor gs:// URIs contain whitespace.
cut -f1,8 "${WORKDIR}/parsed.tsv" | tr '\t' '\n' \
  | xargs -P "${JOBS}" -n2 bash -c 'fetch_one "$1" "$2"' _

GOT="$(find "${STAGE}" -type f \( -name stdout -o -name stderr -o -name rc \) | wc -l | tr -d ' ')"
N="$(wc -l < "${WORKDIR}/parsed.tsv" | tr -d ' ')"
echo "Downloaded ${GOT}/${N} objects." >&2

# Non-zero exit codes, read from the rc files just fetched.
find "${STAGE}" -type f -name rc | while read -r f; do
  v="$(tr -d '[:space:]' < "${f}")"
  [[ "${v}" == "0" || -z "${v}" ]] || printf '%s\t%s\n' "${v}" "${f#${STAGE}/}"
done | sort > "${STAGE}/nonzero_rc.tsv" || true
if [[ -s "${STAGE}/nonzero_rc.tsv" ]]; then
  echo "" >&2
  echo "Non-zero exit codes ($(wc -l < "${STAGE}/nonzero_rc.tsv" | tr -d ' ')):" >&2
  head -20 "${STAGE}/nonzero_rc.tsv" >&2
fi

tar -czf "${OUT}.tar.gz" -C "${WORKDIR}" "${OUT}"
echo "" >&2
echo "Wrote ${OUT}.tar.gz ($(du -h "${OUT}.tar.gz" | cut -f1), ${GOT} objects)" >&2