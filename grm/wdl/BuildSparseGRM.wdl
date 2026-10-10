version 1.0

# Build one SAIGE sparse GRM + PLINK set, standalone, for reuse by every
# association workflow in this repo (FelixPilot, FelixGenome, SaigePilot).
#
# Why this exists: those three workflows each call BuildSaigePlinkAndSparseGRM
# themselves, so "the GRM" is currently whatever a given submission happened to
# build. Running it once here makes it a single versioned artifact with its own
# provenance, and removes the dependence on Cromwell call caching for reuse.
#
# The task below is copied VERBATIM from felix/wdl/FelixGenome.wdl so the GRM
# is identical in kind to the one those workflows produce. Do not edit it here
# without editing it there.
#
# Root the Terra method config on a SET entity (flare_lai_prod_set) so all 22
# chromosomes arrive as one Array[File]:
#   BuildSparseGRM.grm_vcfs = this.flare_lai_prods.gt_vcf
# Confirm the member attribute name in the Terra data tab; it is the pluralized
# member entity type, and this is the repo's first set-rooted config.
#
# NOTE ON MARKER FILTERING: the build script applies none of its own. SAIGE
# filters internally at minMAFforGRM 0.01 and maxMissingRateforGRM 0.15, but
# there is no biallelic-SNV restriction, no long-range-LD/MHC exclusion and no
# LD pruning. That is safe for FLARE anc_vcf input (already biallelic SNVs at
# gnomAD-LAI panel sites) and NOT safe for gt_vcf / ligated_vcf, which are
# unfiltered conversions of the SHAPEIT4 ligated BCF and carry indels and SVs.
# Use anc_vcf here unless the script gains that filtering.
#
# To switch this workflow to gt_vcf / ligated_vcf, build_saige_plink_and_grm.sh
# needs these four steps, the first three applied PER INPUT VCF and before the
# bcftools concat so peak disk stays bounded:
#
#   1. biallelic SNVs      bcftools view -m2 -M2 -v snps
#   2. missingness         bcftools view -i 'F_MISSING <= 0.05'
#   3. long-range LD       bcftools view -t ^chr6:25000000-34000000
#                          (extended MHC, GRCh38. The 8p23.1 and 17q21.31
#                          inversions belong here too, but only with GRCh38
#                          coordinates someone has checked -- do not guess.)
#   4. LD pruning          plink2 --indep-pairwise 50 5 0.2 then --extract,
#                          after --keep and --maf 0.01 --geno 0.05 on the
#                          analysis-sample subset
#
# Then assert the surviving marker count exceeds num_random_markers, or
# createSparseGRM.R fails late with an unhelpful error.
#
# Deliberately NOT in that list: an HWE filter. In a pan-ancestry cohort a
# pooled HWE test rejects markers for ancestry structure rather than for
# genotyping error; it has to be computed within ancestry group or omitted.
#
# Changing the build script is not free: it is a File input, so its hash is
# part of the call-caching key. Editing it in place invalidates every cached
# GRM and every null fit downstream of one. Add a second script and point this
# workflow's build_saige_grm_script at it instead.

workflow BuildSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File build_saige_grm_script
    File make_plink_keep_script
    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 0
  }

  call BuildSaigePlinkAndSparseGRM as MakeGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      build_script = build_saige_grm_script,
      make_plink_keep_script = make_plink_keep_script,
      relatedness_cutoff = relatedness_cutoff,
      num_random_markers = num_random_markers,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      preemptible = preemptible
  }

  output {
    File plink_bed = MakeGRM.plink_bed
    File plink_bim = MakeGRM.plink_bim
    File plink_fam = MakeGRM.plink_fam
    File sparse_grm_mtx = MakeGRM.sparse_grm_mtx
    File sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids
    File paths_tsv = MakeGRM.paths_tsv
  }
}

task BuildSaigePlinkAndSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File build_script
    File make_plink_keep_script
    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    chmod +x "~{build_script}"
    export MAKE_PLINK_KEEP_PY="~{make_plink_keep_script}"
    "~{build_script}" \
      --analysis-samples "~{analysis_samples}" \
      --relatedness-cutoff ~{relatedness_cutoff} \
      --n-threads ~{cpu} \
      --num-random-markers ~{num_random_markers} \
      --out-prefix saige_grm/plink \
      -- \
      ~{sep=" " grm_vcfs}

    cp -L saige_grm/sparseGRM.mtx saige_sparseGRM.mtx
    cp -L saige_grm/sparseGRM.sampleIDs.txt saige_sparseGRM.sampleIDs.txt
    cp saige_grm/plink.bed saige_plink.bed
    cp saige_grm/plink.bim saige_plink.bim
    cp saige_grm/plink.fam saige_plink.fam
    cp saige_grm/saige_sparse_grm.paths.tsv .
  >>>

  output {
    File plink_bed = "saige_plink.bed"
    File plink_bim = "saige_plink.bim"
    File plink_fam = "saige_plink.fam"
    File sparse_grm_mtx = "saige_sparseGRM.mtx"
    File sparse_grm_sample_ids = "saige_sparseGRM.sampleIDs.txt"
    File paths_tsv = "saige_sparse_grm.paths.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
