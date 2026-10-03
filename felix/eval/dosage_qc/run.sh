#!/usr/bin/env bash
# Runs inside felix-pilot:0.1.0 with the spike dir mounted at /w.
set -euo pipefail
cd /w
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

for f in joint snv_admix rep_s1 rep_s2 rep_shift; do
  bcftools view -Oz -o $f.vcf.gz $f.vcf
  bcftools index -f -c $f.vcf.gz
  bcftools index -f -t $f.vcf.gz
done
plink2 --vcf grm.vcf --make-bed --out grm --silent --allow-extra-chr

mkdir -p out
for ph in y_bin10 y_bin2 y_q; do
  tt=binary; [[ $ph == y_q ]] && tt=quantitative
  if [[ ! -f out/null_$ph.rda ]]; then
    Rscript /usr/local/bin/step1_fitNULLGLMM.R \
      --plinkFile=grm --phenoFile=pheno.tsv --phenoCol=$ph \
      --covarColList=x1,q --sampleIDColinphenoFile=ID --traitType=$tt \
      --outputPrefix=out/null_$ph --nThreads=2 --LOCO=FALSE \
      --IsOverwriteVarianceRatioFile=TRUE --isCateVarianceRatio=FALSE \
      > out/step1_$ph.log 2>&1
  fi
done

cut -f1 pheno.tsv | tail -n +2 > keep.txt
felixla --phase-vcf joint.vcf.gz --flare-vcf joint.vcf.gz --n-ancestries 2 \
  --keep keep.txt --make-felixla --out out/packed > out/felixla.log 2>&1

common=(--chrom=chr1 --is_admixed=TRUE --number_of_ancestry=2
        --pvalcutoff_of_haplotype=0.05 --minMAC=1 --nThreads=1
        --is_Firth_beta=TRUE --is_output_moreDetails=TRUE --LOCO=FALSE)
for ph in y_bin10 y_bin2 y_q; do
  m=(--GMMATmodelFile=out/null_$ph.rda --varianceRatioFile=out/null_$ph.varianceRatio.txt)
  Rscript /usr/local/bin/step2_SPAtests.R --FELIXlaPrefix=out/packed "${m[@]}" "${common[@]}" \
    --SAIGEOutputFile=out/felixla_$ph.txt > out/step2_felixla_$ph.log 2>&1 || echo "FAIL felixla $ph"
  for v in snv_admix rep_s1 rep_s2 rep_shift; do
    Rscript /usr/local/bin/step2_SPAtests.R --vcfFile=$v.vcf.gz --vcfFileIndex=$v.vcf.gz.csi \
      --vcfField=DS "${m[@]}" "${common[@]}" \
      --SAIGEOutputFile=out/${v}_$ph.txt > out/step2_${v}_$ph.log 2>&1 || echo "FAIL $v $ph"
  done
done
ls -la out | grep -v "\.log$"
