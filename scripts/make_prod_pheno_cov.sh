#!/usr/bin/env bash
# Build and upload the four restricted pheno_cov files the felix.prod.* run
# table expects. Run on the Terra VM, where the real pheno_cov.tsv can be
# fetched; it cannot be rebuilt from the source tables without risking
# covariate values that differ from the 947 phenotypes using the shared file.
#
# PhecodeX labels all four "Both", but their cases are (almost) all female, so
# the sex term is separated and step 1 diverges -- which is what produced the
# gamma_q_inv overflow in step 2. Setting the excluded sex to NA in that one
# phenotype's column makes SAIGE drop those participants for that phenotype
# only. analysis_samples is deliberately NOT touched: it feeds PackFelixla and
# MakeGRM, and changing it would re-pack 22 chromosomes.
set -euo pipefail

: "${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET}"
SHARED="${WORKSPACE_BUCKET}/felix/phenotypes/wgs/pheno_cov.tsv"
DEST="${WORKSPACE_BUCKET}/felix/phenotypes/prod"
COV="${COV:-covariates.v8.csv.gz}"
PHENOS=(CA_142_2 CA_105 GU_614_55 CA_142_21)

[[ -f scripts/make_restricted_pheno_cov.py ]] || {
  echo "run from the repo root after 00_sync_repo" >&2; exit 1; }
[[ -f "$COV" ]] || { echo "missing $COV (set COV=/path/to/covariates.v8.csv.gz)" >&2; exit 1; }

[[ -f pheno_cov.tsv ]] || { echo "fetching $SHARED"; gsutil cp "$SHARED" pheno_cov.tsv; }

for p in "${PHENOS[@]}"; do
  echo "=== $p ==="
  # Output name carries the BARE phenotype, not the felix.prod. row id: that is
  # what the four entity rows reference. A mismatch fails at localization
  # inside the workflow, not here.
  python3 scripts/make_restricted_pheno_cov.py \
    --pheno-cov pheno_cov.tsv --covariates "$COV" \
    --phenotype "$p" --restrict-to Female \
    --out "pheno_cov.$p.tsv"
done

echo
echo "=== row counts (must all equal the shared file) ==="
shared_rows=$(wc -l < pheno_cov.tsv | tr -d ' ')
echo "pheno_cov.tsv $shared_rows"
for p in "${PHENOS[@]}"; do
  n=$(wc -l < "pheno_cov.$p.tsv" | tr -d ' ')
  echo "pheno_cov.$p.tsv $n"
  [[ "$n" == "$shared_rows" ]] || { echo "  ROW COUNT CHANGED -- do not upload" >&2; exit 1; }
done

echo
echo "=== uploading to $DEST/ ==="
for p in "${PHENOS[@]}"; do gsutil cp "pheno_cov.$p.tsv" "$DEST/"; done
gsutil ls "$DEST/pheno_cov.*.tsv"
