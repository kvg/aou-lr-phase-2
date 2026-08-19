#!/usr/bin/env bash
# Build a PLINK bed dataset from GRM VCFs and a SAIGE sparse GRM.
#
# Inputs (env or flags):
#   --analysis-samples  sample ID list (VCF / analysis order)
#   --grm-vcfs          space-separated VCF paths (or pass as trailing args after --)
#   --relatedness-cutoff  SAIGE relatednessCutoff (default 0.05)
#   --n-threads
#   --out-prefix        output prefix directory/prefix (default: saige_grm/plink)
#
# Writes:
#   ${prefix}.bed/.bim/.fam
#   sparse GRM mtx + sampleIDs from createSparseGRM.R
#   saige_sparse_grm.paths.tsv documenting the produced file names

set -euo pipefail

ANALYSIS_SAMPLES=""
RELATEDNESS_CUTOFF="0.05"
N_THREADS="8"
OUT_PREFIX="saige_grm/plink"
NUM_RANDOM_MARKERS="2000"
CREATE_SPARSE_GRM_R="${CREATE_SPARSE_GRM_R:-/opt/SAIGE/extdata/createSparseGRM.R}"
VCFS=()

usage() {
  cat <<EOF
Usage: $0 --analysis-samples FILE --out-prefix PREFIX [options] -- vcf1 [vcf2 ...]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-samples) ANALYSIS_SAMPLES="$2"; shift 2 ;;
    --relatedness-cutoff) RELATEDNESS_CUTOFF="$2"; shift 2 ;;
    --n-threads) N_THREADS="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --num-random-markers) NUM_RANDOM_MARKERS="$2"; shift 2 ;;
    --create-sparse-grm-r) CREATE_SPARSE_GRM_R="$2"; shift 2 ;;
    --) shift; VCFS+=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *)
      if [[ -f "$1" ]]; then
        VCFS+=("$1"); shift
      else
        echo "Unknown arg: $1" >&2; usage >&2; exit 1
      fi
      ;;
  esac
done

if [[ -z "${ANALYSIS_SAMPLES}" || ${#VCFS[@]} -eq 0 ]]; then
  echo "Required: --analysis-samples and at least one VCF" >&2
  usage >&2
  exit 1
fi
if [[ ! -x "${CREATE_SPARSE_GRM_R}" && ! -f "${CREATE_SPARSE_GRM_R}" ]]; then
  echo "createSparseGRM.R not found: ${CREATE_SPARSE_GRM_R}" >&2
  exit 1
fi

OUT_DIR="$(dirname "${OUT_PREFIX}")"
mkdir -p "${OUT_DIR}" grm_build
VCF_LIST=grm_build/vcf_list.txt
: > "${VCF_LIST}"
for vcf in "${VCFS[@]}"; do
  echo "${vcf}" >> "${VCF_LIST}"
done

N_VCF=$(wc -l < "${VCF_LIST}" | tr -d ' ')
if [[ "${N_VCF}" -eq 1 ]]; then
  INPUT_VCF=$(head -n1 "${VCF_LIST}")
else
  bcftools concat -f "${VCF_LIST}" -Oz -o grm_build/concat.vcf.gz --threads "${N_THREADS}"
  bcftools index -t grm_build/concat.vcf.gz
  INPUT_VCF=grm_build/concat.vcf.gz
fi

plink2 --vcf "${INPUT_VCF}" \
  --double-id \
  --make-bed \
  --set-all-var-ids "@:#:\$r:\$a" \
  --out grm_build/all \
  --threads "${N_THREADS}"

KEEP_PY="${MAKE_PLINK_KEEP_PY:-$(dirname "${BASH_SOURCE[0]}")/make_plink_keep.py}"
if [[ ! -f "${KEEP_PY}" ]]; then
  KEEP_PY="/opt/tractor_mix_scripts/make_plink_keep.py"
fi
python3 "${KEEP_PY}" \
  --analysis-samples "${ANALYSIS_SAMPLES}" \
  --fam grm_build/all.fam \
  --out-keep grm_build/keep_iid.txt \
  --out-samples grm_build/analysis_samples.intersect.txt

plink2 --bfile grm_build/all \
  --keep grm_build/keep_iid.txt \
  --make-bed \
  --out "${OUT_PREFIX}" \
  --threads "${N_THREADS}"

# SAIGE createSparseGRM expects a PLINK prefix without extension
Rscript "${CREATE_SPARSE_GRM_R}" \
  --plinkFile="${OUT_PREFIX}" \
  --nThreads="${N_THREADS}" \
  --outputPrefix="${OUT_DIR}/sparseGRM" \
  --numRandomMarkerforSparseKin="${NUM_RANDOM_MARKERS}" \
  --relatednessCutoff="${RELATEDNESS_CUTOFF}"

# Locate outputs (SAIGE embeds cutoff + marker count in the filename)
SPARSE_MTX=$(ls -1 "${OUT_DIR}"/sparseGRM*.sparseGRM.mtx | head -n1)
SPARSE_IDS=$(ls -1 "${OUT_DIR}"/sparseGRM*.sparseGRM.mtx.sampleIDs.txt | head -n1)
if [[ -z "${SPARSE_MTX}" || -z "${SPARSE_IDS}" ]]; then
  echo "Failed to locate createSparseGRM outputs under ${OUT_DIR}" >&2
  ls -la "${OUT_DIR}" >&2 || true
  exit 1
fi

# Stable symlinks for WDL outputs
ln -sf "$(basename "${SPARSE_MTX}")" "${OUT_DIR}/sparseGRM.mtx"
ln -sf "$(basename "${SPARSE_IDS}")" "${OUT_DIR}/sparseGRM.sampleIDs.txt"

cat > "${OUT_DIR}/saige_sparse_grm.paths.tsv" <<EOF
key	path
plink_prefix	${OUT_PREFIX}
sparse_grm_mtx	${SPARSE_MTX}
sparse_grm_sample_ids	${SPARSE_IDS}
relatedness_cutoff	${RELATEDNESS_CUTOFF}
num_random_markers	${NUM_RANDOM_MARKERS}
EOF

echo "Wrote PLINK prefix ${OUT_PREFIX} and sparse GRM:"
echo "  ${SPARSE_MTX}"
echo "  ${SPARSE_IDS}"
