#!/bin/bash
# Locityper container worker (this image has locityper + bash + GNU parallel, not python3).
# Host already downloaded reads/reference/counts/BED/db onto /work.
set -euo pipefail
cd /work
if [ -f /work/task.env ]; then
  # shellcheck disable=SC1091
  . /work/task.env
fi
: "${SAMPLE_ID:?SAMPLE_ID empty}"
: "${LOCITYPER_N_CPU:=2}"
: "${TECHNOLOGY:=illumina}"
export SHELL=/bin/bash

exec >>/work/worker.log 2>&1
echo "genotype start sample=${SAMPLE_ID} n_cpu=${LOCITYPER_N_CPU} tech=${TECHNOLOGY}"

(
  printf 'timestamp\trss_kb\n' > /work/resource_stats.tsv
  while true; do
    rss=$(awk '/^VmRSS:/{s+=$2} END {print s+0}' /proc/[0-9]*/status 2>/dev/null || echo 0)
    printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${rss}" >> /work/resource_stats.tsv
    sleep 5
  done
) &
MON_PID=$!
trap 'kill "${MON_PID}" 2>/dev/null || true' EXIT

for f in reads.cram reads.cram.crai reference.fa reference.fa.fai counts.jf loci.bed vcf_db.tar.gz; do
  if [ ! -s "${f}" ]; then
    echo "missing input: /work/${f}" >&2
    exit 1
  fi
done

ln -sfn reads.cram subset.cram
ln -sfn reads.cram.crai subset.cram.crai

if ! command -v locityper >/dev/null 2>&1; then
  echo "locityper is not on PATH" >&2
  exit 1
fi

echo "+ locityper preproc ..."
locityper preproc \
  -a subset.cram \
  -r reference.fa \
  -j counts.jf \
  -@ "${LOCITYPER_N_CPU}" \
  --technology "${TECHNOLOGY}" \
  -o locityper_preproc

tar -xzf vcf_db.tar.gz
if [ ! -d vcf_db ]; then
  echo "vcf_db.tar.gz did not contain vcf_db/" >&2
  exit 1
fi

mkdir -p out_dir/loci

process_single_locus() {
  line="$1"
  case "${line}" in
    ""|\#*) return 0 ;;
  esac
  locus_name=$(printf '%s\n' "${line}" | cut -f4)
  if [ -z "${locus_name}" ]; then
    echo "BED line missing column 4: ${line}" >&2
    return 1
  fi
  mkdir -p "out_dir/loci/${locus_name}"
  locityper genotype \
    -a subset.cram \
    -r reference.fa \
    -d vcf_db \
    -p locityper_preproc \
    --subset-loci "${locus_name}" \
    -o out_dir >/dev/null 2>&1
}
export -f process_single_locus

if [ ! -x /usr/bin/parallel ]; then
  echo "GNU parallel missing at /usr/bin/parallel" >&2
  command -v parallel || true
  exit 1
fi

grep -v '^#' loci.bed | grep -ve '^[[:space:]]*$' > /work/loci.bed.nz
echo "loci: $(wc -l < /work/loci.bed.nz)"
/usr/bin/parallel --line-buffer -j "${LOCITYPER_N_CPU}" process_single_locus {} < /work/loci.bed.nz

find out_dir -type f -name '*.bam' -exec rm -f {} \;
tar -czf "${SAMPLE_ID}.locityper.tar.gz" out_dir
ls -lh "${SAMPLE_ID}.locityper.tar.gz"
echo "genotype container done"
