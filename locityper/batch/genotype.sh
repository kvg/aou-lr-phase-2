#!/usr/bin/env bash
# Cloud Batch / dsub: locityper preproc + genotype (all BED loci) + summarize.
# No BED scatter — GNU parallel on this VM, same as the WDL shard internals.
set -euo pipefail
date

: "${SAMPLE_ID:?}"
: "${MINICRAM:?}"
: "${MINICRAI:?}"
: "${REF_FA:?}"
: "${REF_FAI:?}"
: "${COUNTS_JF:?}"
: "${BED:?}"
: "${DB_TAR:?}"
: "${SUMMARY_CSV:?}"
: "${RESULTS_TAR:?}"

N_CPU="${LOCITYPER_N_CPU:-2}"
TECHNOLOGY="${TECHNOLOGY:-illumina}"

ln -sf "${REF_FA}" reference.fa
ln -sf "${REF_FAI}" reference.fa.fai
ln -sf "${MINICRAM}" subset.cram
ln -sf "${MINICRAI}" subset.cram.crai

locityper preproc -a subset.cram \
  -r reference.fa \
  -j "${COUNTS_JF}" \
  -@ "${N_CPU}" \
  --technology "${TECHNOLOGY}" \
  -o locityper_preproc

tar -xzf "${DB_TAR}"
mkdir -p out_dir/loci

process_single_locus() {
  line="$1"
  locus_name=$(echo "$line" | cut -f4)
  if [[ -z "${locus_name}" ]]; then
    echo "BED line missing column 4 (locus name): ${line}" >&2
    return 1
  fi
  mkdir -p "out_dir/loci/${locus_name}"
  locityper genotype -a subset.cram \
    -r reference.fa \
    -d vcf_db \
    -p locityper_preproc \
    --subset-loci "${locus_name}" \
    -o out_dir >/dev/null 2>&1
}
export -f process_single_locus

cat "${BED}" | /usr/bin/parallel --line-buffer -j "${N_CPU}" process_single_locus {}
find out_dir -type f -name "*.bam" -exec rm -f {} \;

tar -czf "${SAMPLE_ID}.locityper.tar.gz" out_dir

python3 - <<PY
import gzip, json, math
from pathlib import Path

sample_id = "${SAMPLE_ID}"
loci_root = Path("out_dir") / "loci"
out_path = Path("gts.filtered.csv")
rows = ["sample\tlocus\tgenotype\tquality\ttotal_reads\tunexpl_reads\tweight_dist\twarnings"]
if loci_root.is_dir():
    for entry in sorted(loci_root.iterdir(), key=lambda p: p.name):
        if not entry.is_dir():
            continue
        json_path = entry / "res.json.gz"
        if not json_path.is_file():
            continue
        with gzip.open(json_path, "rt") as fh:
            res = json.load(fh)
        line = f"{sample_id}\t{entry.name}\t"
        if "genotype" not in res:
            rows.append(line + "*")
            continue
        gt = res["genotype"]
        qual = math.floor(10 * float(res["quality"])) * 0.1
        total_reads = res.get("total_reads", "")
        unexpl_reads = res.get("unexpl_reads", "")
        weight_dist = res.get("weight_dist")
        weight_s = "" if weight_dist is None else f"{float(weight_dist):.5f}"
        warnings = ";".join(res.get("warnings", [])) or "*"
        rows.append(
            line + f"{gt}\t{qual:.1f}\t{total_reads}\t{unexpl_reads}\t{weight_s}\t{warnings}"
        )
out_path.write_text("\n".join(rows) + "\n")
print(f"wrote {len(rows) - 1} loci to {out_path}", flush=True)
PY

mkdir -p "$(dirname "${SUMMARY_CSV}")" "$(dirname "${RESULTS_TAR}")"
cp -f gts.filtered.csv "${SUMMARY_CSV}"
cp -f "${SAMPLE_ID}.locityper.tar.gz" "${RESULTS_TAR}"
ls -lh "${SUMMARY_CSV}" "${RESULTS_TAR}"
date
