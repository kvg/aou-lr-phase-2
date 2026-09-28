#!/usr/bin/env bash
# Cloud Batch / dsub: ExpansionHunter on a local minicram.
set -euo pipefail
date

: "${SAMPLE_ID:?}"
: "${MINICRAM:?}"
: "${MINICRAI:?}"
: "${REF_FA:?}"
: "${REF_FAI:?}"
: "${CATALOG:?}"
: "${EH_JSON:?}"
: "${EH_VCF:?}"
SEX="${SEX:-female}"

command -v ExpansionHunter
ExpansionHunter --version 2>&1 || true
ExpansionHunter --help > eh.help.txt 2>&1 || true

sex_lc=$(printf '%s' "${SEX}" | tr '[:upper:]' '[:lower:]' | tr -d ' ')
case "${sex_lc}" in
  male|m|1)   sex_flag="male"   ;;
  female|f|2) sex_flag="female" ;;
  *)
    echo "WARN: unrecognized sex '${SEX}'; defaulting to female" >&2
    sex_flag="female"
    ;;
esac
echo "sex=${sex_flag}"

ln -sf "${REF_FA}" reference.fa
ln -sf "${REF_FAI}" reference.fa.fai
ln -sf "${MINICRAM}" reads.cram
ln -sf "${MINICRAI}" reads.cram.crai

if grep -q -- '--variant-catalog' eh.help.txt; then
  CAT_FLAG="--variant-catalog"
else
  CAT_FLAG="--catalog"
fi
if grep -q -- 'optimized-streaming' eh.help.txt; then
  MODE="optimized-streaming"
else
  MODE="seeking"
fi
echo "catalog flag=${CAT_FLAG} analysis-mode=${MODE}"

extra=()
grep -q -- '--dont-output-consensus-sequences' eh.help.txt \
  && extra+=(--dont-output-consensus-sequences)
grep -q -- '--dont-output-quality-metrics' eh.help.txt \
  && extra+=(--dont-output-quality-metrics)
grep -q -- '--disable-all-plots' eh.help.txt \
  && extra+=(--disable-all-plots)

prefix="${SAMPLE_ID}.EH"
ExpansionHunter \
  --reads reads.cram \
  --reads-index reads.cram.crai \
  --reference reference.fa \
  "${CAT_FLAG}" "${CATALOG}" \
  --output-prefix "${prefix}" \
  --sex "${sex_flag}" \
  --analysis-mode "${MODE}" \
  "${extra[@]}"

if [[ ! -s "${prefix}.json" || ! -s "${prefix}.vcf" ]]; then
  echo "ExpansionHunter missing ${prefix}.json / ${prefix}.vcf" >&2
  ls -la
  exit 1
fi
python3 -c "import json; json.load(open('${prefix}.json'))"

mkdir -p "$(dirname "${EH_JSON}")" "$(dirname "${EH_VCF}")"
cp -f "${prefix}.json" "${EH_JSON}"
cp -f "${prefix}.vcf" "${EH_VCF}"
ls -lh "${EH_JSON}" "${EH_VCF}"
date
