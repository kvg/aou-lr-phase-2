version 1.0

# Self-contained SAIGE chr22 calibration pilot for Terra Methods Repository.
# Matched interface to TractorMixPilot: shared pheno_cov / analysis_samples /
# selected_phenotypes / covariate_columns (limited or recommended full).

task BuildSaigePlinkAndSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File build_script
    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 2
  }

  command <<<
    set -euo pipefail
    chmod +x "~{build_script}"
    "~{build_script}" \
      --analysis-samples "~{analysis_samples}" \
      --relatedness-cutoff ~{relatedness_cutoff} \
      --n-threads ~{cpu} \
      --num-random-markers ~{num_random_markers} \
      --out-prefix saige_grm/plink \
      -- \
      ~{sep=" " grm_vcfs}

    # Materialize symlink targets as real files for Cromwell localization
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

task FitSaigeNull {
  input {
    File pheno_cov
    String phenotype
    File covariate_columns
    File analysis_samples
    File plink_bed
    File plink_bim
    File plink_fam
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File fit_null_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 8
    Int memory_gb = 16
    Int disk_gb = 100
    Int preemptible = 2
  }

  command <<<
    set -euo pipefail
    mkdir -p null
    # Reconstruct PLINK prefix beside localized bed/bim/fam
    ln -sf "~{plink_bed}" null/plink.bed
    ln -sf "~{plink_bim}" null/plink.bim
    ln -sf "~{plink_fam}" null/plink.fam

    Rscript "~{fit_null_script}" \
      --pheno-cov "~{pheno_cov}" \
      --phenotype "~{phenotype}" \
      --covariates "~{covariate_columns}" \
      --analysis-samples "~{analysis_samples}" \
      --plink-prefix null/plink \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --n-threads ~{cpu} \
      --out-prefix "null/~{phenotype}"

    # Collect SAIGE null artifacts under stable names
    cp "null/~{phenotype}.rda" "~{phenotype}.null.rda" || \
      cp $(ls -1 null/~{phenotype}*.rda | head -n1) "~{phenotype}.null.rda"
    cp "null/~{phenotype}.varianceRatio.txt" "~{phenotype}.varianceRatio.txt" || \
      cp $(ls -1 null/~{phenotype}*.varianceRatio.txt | head -n1) "~{phenotype}.varianceRatio.txt"
    cp "null/~{phenotype}.pheno.tsv" "~{phenotype}.pheno.tsv"
    cp "null/~{phenotype}.samples.txt" "~{phenotype}.samples.txt"
    cp "null/~{phenotype}.null_meta.tsv" "~{phenotype}.null_meta.tsv"
  >>>

  output {
    File null_rda = "~{phenotype}.null.rda"
    File variance_ratio = "~{phenotype}.varianceRatio.txt"
    File pheno_used = "~{phenotype}.pheno.tsv"
    File samples_used = "~{phenotype}.samples.txt"
    File null_meta = "~{phenotype}.null_meta.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task RunSaigeStep2 {
  input {
    File flare_vcf
    String chrom = "chr22"
    String phenotype
    File null_rda
    File variance_ratio
    File samples_used
    File run_step2_script
    Int min_mac = 20
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 8
    Int memory_gb = 16
    Int disk_gb = 200
    Int preemptible = 2
  }

  command <<<
    set -euo pipefail
    mkdir -p step2
    # Reconstruct null prefix expected by run_saige_step2.R
    ln -sf "~{null_rda}" "step2/~{phenotype}.rda"
    ln -sf "~{variance_ratio}" "step2/~{phenotype}.varianceRatio.txt"

    # Localize / index VCF
    VCF="~{flare_vcf}"
    if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
      bcftools index -t "${VCF}" || bcftools index -c "${VCF}"
    fi

    Rscript "~{run_step2_script}" \
      --vcf "${VCF}" \
      --chrom "~{chrom}" \
      --null-prefix "step2/~{phenotype}" \
      --sample-file "~{samples_used}" \
      --min-mac ~{min_mac} \
      --n-threads ~{cpu} \
      --out-tsv "~{phenotype}.saige.tsv" \
      --out-raw "~{phenotype}.saige.raw.txt"
  >>>

  output {
    File results_tsv = "~{phenotype}.saige.tsv"
    File results_raw = "~{phenotype}.saige.raw.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

workflow SaigePilot {
  input {
    File flare_vcf
    String chrom = "chr22"

    # Shared cohort / table from notebooks/tractor_01_prepare_inputs.ipynb.
    # `pheno_cov` contains phenotypes and both covariate matrices; select one
    # with covariate_columns_limited.txt or covariate_columns_full.txt.
    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns

    Array[File] grm_vcfs

    File build_saige_grm_script
    File fit_saige_null_script
    File run_saige_step2_script

    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    Int min_mac = 20

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)

  call BuildSaigePlinkAndSparseGRM as BuildGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      build_script = build_saige_grm_script,
      relatedness_cutoff = relatedness_cutoff,
      num_random_markers = num_random_markers,
      docker = docker
  }

  scatter (pheno in phenotypes) {
    call FitSaigeNull as Null {
      input:
        pheno_cov = pheno_cov,
        phenotype = pheno,
        covariate_columns = covariate_columns,
        analysis_samples = analysis_samples,
        plink_bed = BuildGRM.plink_bed,
        plink_bim = BuildGRM.plink_bim,
        plink_fam = BuildGRM.plink_fam,
        sparse_grm_mtx = BuildGRM.sparse_grm_mtx,
        sparse_grm_sample_ids = BuildGRM.sparse_grm_sample_ids,
        fit_null_script = fit_saige_null_script,
        docker = docker
    }

    call RunSaigeStep2 as Score {
      input:
        flare_vcf = flare_vcf,
        chrom = chrom,
        phenotype = pheno,
        null_rda = Null.null_rda,
        variance_ratio = Null.variance_ratio,
        samples_used = Null.samples_used,
        run_step2_script = run_saige_step2_script,
        min_mac = min_mac,
        docker = docker
    }
  }

  output {
    File plink_bed = BuildGRM.plink_bed
    File plink_bim = BuildGRM.plink_bim
    File plink_fam = BuildGRM.plink_fam
    File sparse_grm_mtx = BuildGRM.sparse_grm_mtx
    File sparse_grm_sample_ids = BuildGRM.sparse_grm_sample_ids
    Array[File] null_rdas = Null.null_rda
    Array[File] null_metas = Null.null_meta
    Array[File] results_tsvs = Score.results_tsv
    Array[File] results_raw = Score.results_raw
  }

  meta {
    description: "SAIGE chr22 pilot matched to Tractor-Mix: sparse GRM from chr1+chr22, per-phenotype logistic null + SPA tests."
    allowNestedInputs: true
  }
}
