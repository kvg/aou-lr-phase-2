#!/usr/bin/env bash
# Runs inside felix-dosageqc-test:local with the spike dir at /w.
set -uo pipefail
cd /w
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p val
common=(--chrom=chr1 --is_admixed=TRUE --number_of_ancestry=2 --pvalcutoff_of_haplotype=0.05
        --minMAC=1 --nThreads=1 --is_Firth_beta=TRUE --is_output_moreDetails=TRUE --LOCO=FALSE)
m() { echo --GMMATmodelFile=out/null_$1.rda --varianceRatioFile=out/null_$1.varianceRatio.txt; }

# V1: switch off -> must match stock FELIX outputs in out/.
unset FELIX_DOSAGE_QC
Rscript /usr/local/bin/step2_SPAtests.R --FELIXlaPrefix=out/packed $(m y_bin10) "${common[@]}" \
  --SAIGEOutputFile=val/off_felixla_y_bin10.txt > val/off_felixla.log 2>&1 || echo FAIL v1a
for v in snv_admix rep_s1; do
  Rscript /usr/local/bin/step2_SPAtests.R --vcfFile=$v.vcf.gz --vcfFileIndex=$v.vcf.gz.csi --vcfField=DS \
    $(m y_bin10) "${common[@]}" --SAIGEOutputFile=val/off_${v}_y_bin10.txt > val/off_$v.log 2>&1 || echo FAIL v1 $v
done

# V2/V3: switch on.
export FELIX_DOSAGE_QC=carrier
Rscript /usr/local/bin/step2_SPAtests.R --vcfFile=snv_admix.vcf.gz --vcfFileIndex=snv_admix.vcf.gz.csi --vcfField=DS \
  $(m y_bin10) "${common[@]}" --SAIGEOutputFile=val/on_snv_admix_y_bin10.txt > val/on_snv.log 2>&1 || echo FAIL v2
for ph in y_bin2 y_bin10 y_q; do
  Rscript /usr/local/bin/step2_SPAtests.R --vcfFile=rep_null.vcf.gz --vcfFileIndex=rep_null.vcf.gz.csi --vcfField=DS \
    $(m $ph) "${common[@]}" --minMAC=20 --SAIGEOutputFile=val/on_rep_null_$ph.txt > val/on_rep_null_$ph.log 2>&1 || echo FAIL v3 $ph
done
unset FELIX_DOSAGE_QC

grep -h "FELIX_DOSAGE_QC=carrier" val/*.log | sort | uniq -c
echo done
