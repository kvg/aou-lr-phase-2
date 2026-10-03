#!/usr/bin/env bash
set -euo pipefail
cd /w
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
bcftools view -Oz -o rep_null.vcf.gz rep_null.vcf && bcftools index -f -c rep_null.vcf.gz
common=(--chrom=chr1 --is_admixed=TRUE --number_of_ancestry=2 --pvalcutoff_of_haplotype=0.05
        --minMAC=1 --nThreads=1 --is_Firth_beta=TRUE --is_output_moreDetails=TRUE --LOCO=FALSE)
for ph in y_bin2 y_bin10 y_q; do
  Rscript /usr/local/bin/step2_SPAtests.R --vcfFile=rep_null.vcf.gz --vcfFileIndex=rep_null.vcf.gz.csi \
    --vcfField=DS --GMMATmodelFile=out/null_$ph.rda --varianceRatioFile=out/null_$ph.varianceRatio.txt \
    "${common[@]}" --SAIGEOutputFile=out/rep_null_$ph.txt > out/step2_rep_null_$ph.log 2>&1 || echo "FAIL $ph"
done
echo done
