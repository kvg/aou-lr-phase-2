#!/usr/bin/env bash
set -uo pipefail
cd /w
mkdir -p val
# The wrapper only checks that the sparse GRM files exist; the stub never reads them.
touch val/sparse.mtx val/sparse.ids
touch val/stub.rda val/stub.varianceRatio.txt
export FELIX_DOSAGE_QC=leaked
Rscript /repo/felix/scripts/run_felix_step2.R --vcf-file rep_null.vcf.gz --dosage-qc carrier --chrom chr1 \
  --null-prefix val/stub --sparse-grm val/sparse.mtx --sparse-grm-ids val/sparse.ids \
  --min-mac 20 --n-ancestries 2 --step2-r /repo/felix/eval/dosage_qc/stub_step2.R --out-tsv val/stub_out.tsv > val/stub_run1.log 2>&1 && echo "run1 ok"
cp val/stub_args.txt val/stub_args_vcf.txt
Rscript /repo/felix/scripts/run_felix_step2.R --felixla-prefix out/packed --chrom chr1 \
  --null-prefix val/stub --sparse-grm val/sparse.mtx --sparse-grm-ids val/sparse.ids \
  --n-ancestries 2 --step2-r /repo/felix/eval/dosage_qc/stub_step2.R --out-tsv val/stub_out2.tsv > val/stub_run2.log 2>&1 && echo "run2 ok"
cp val/stub_args.txt val/stub_args_felixla.txt
Rscript /repo/felix/scripts/run_felix_step2.R --felixla-prefix out/packed --dosage-qc carrier --chrom chr1 --null-prefix val/stub \
  --sparse-grm val/sparse.mtx --sparse-grm-ids val/sparse.ids --step2-r /repo/felix/eval/dosage_qc/stub_step2.R --out-tsv val/x.tsv 2>&1 | grep -m1 -o "applies to --vcf-file input only"
Rscript /repo/felix/scripts/run_felix_step2.R --felixla-prefix out/packed --vcf-file rep_null.vcf.gz --chrom chr1 --null-prefix val/stub \
  --sparse-grm val/sparse.mtx --sparse-grm-ids val/sparse.ids --step2-r /repo/felix/eval/dosage_qc/stub_step2.R --out-tsv val/x.tsv 2>&1 | grep -m1 -o "exactly one of"
echo "--- vcf mode args"; cat val/stub_args_vcf.txt
echo "--- felixla mode args"; cat val/stub_args_felixla.txt
