version 1.0

# Checks per-locus repeat dosage from the integrated callset against TRGT on
# one chromosome for a participant subset: exact / within-one-unit agreement
# and integrated-sum overcounts (split or double-counted records), stratified
# by how many integrated records hit the locus. See
# scripts/compare_trgt_locus_dosage.py.
#
# trgt_table: participant<TAB>TRGT VCF URI, the subset to compare.
# vcfs / vcf_roles: the same callsets and roles as RepeatLocusDosage
# (Phase 2: GLnexus shard for `chrom` as "small", v3_main.bcf as "sv").

workflow TrgtLocusConcordance {
  input {
    String chrom = "chr20"
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] vcf_roles
    File catalog_bed
    File trgt_table
    Int flank = 10
    Int split_bp = 50
    String apply_filters = "PASS,."
    Boolean ignore_phase = true

    File aggregate_repeat_loci_py
    File compare_trgt_locus_dosage_py

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int memory_gb = 16
    Int preemptible = 1
  }

  Array[Array[String]] columns = transpose(read_tsv(trgt_table))

  call Concordance {
    input:
      chrom = chrom,
      vcfs = vcfs,
      vcf_indexes = vcf_indexes,
      vcf_roles = vcf_roles,
      split_bp = split_bp,
      apply_filters = apply_filters,
      ignore_phase = ignore_phase,
      catalog_bed = catalog_bed,
      samples = columns[0],
      trgt_vcfs = columns[1],
      flank = flank,
      aggregate_repeat_loci_py = aggregate_repeat_loci_py,
      compare_trgt_locus_dosage_py = compare_trgt_locus_dosage_py,
      docker = docker,
      memory_gb = memory_gb,
      preemptible = preemptible
  }

  output {
    File loci = Concordance.loci
    File summary = Concordance.summary
    File aggregate_summary = Concordance.aggregate_summary
  }

  meta {
    description: "Integrated-callset locus dosage vs TRGT on one chromosome (double-count check for Figure 2C)."
    allowNestedInputs: true
  }
}

task Concordance {
  input {
    String chrom
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] vcf_roles
    Int split_bp
    String apply_filters
    Boolean ignore_phase
    File catalog_bed
    Array[String] samples
    Array[File] trgt_vcfs
    Int flank
    File aggregate_repeat_loci_py
    File compare_trgt_locus_dosage_py
    String docker
    Int memory_gb
    Int preemptible
  }

  Int disk_gb = ceil((size(vcfs, "GB") + size(trgt_vcfs, "GB")) * 1.2) + 30

  command <<<
    set -euo pipefail
    mkdir -p in
    paste "~{write_lines(vcfs)}" "~{write_lines(vcf_indexes)}" "~{write_lines(vcf_roles)}" > pairs.tsv
    VCF_ARGS=()
    n=0
    while IFS=$'\t' read -r v x role; do
      n=$((n + 1))
      name="in/${n}.$(basename "${v}")"
      ln -s "${v}" "${name}"
      case "${x}" in *.csi) ln -s "${x}" "${name}.csi" ;; *) ln -s "${x}" "${name}.tbi" ;; esac
      case "${role}" in
        all) VCF_ARGS+=(--vcf "${name}") ;;
        small) VCF_ARGS+=(--small-vcf "${name}") ;;
        sv) VCF_ARGS+=(--sv-vcf "${name}") ;;
        *) echo "unknown VCF role ${role}" >&2; exit 1 ;;
      esac
    done < pairs.tsv
    cp "~{write_lines(samples)}" samples.txt
    paste samples.txt "~{write_lines(trgt_vcfs)}" > trgt.tsv

    python3 "~{aggregate_repeat_loci_py}" \
      "${VCF_ARGS[@]}" \
      --catalog-bed "~{catalog_bed}" \
      --chrom "~{chrom}" \
      --samples samples.txt \
      --flank ~{flank} \
      --split-bp ~{split_bp} \
      --apply-filters "~{apply_filters}" \
      ~{true="--ignore-phase" false="" ignore_phase} \
      --out-alleles alleles.tsv.gz \
      --out-vcf loci.vcf \
      --out-summary aggregate.summary.json

    python3 "~{compare_trgt_locus_dosage_py}" \
      --locus-vcf loci.vcf \
      --catalog-bed "~{catalog_bed}" \
      --trgt-tsv trgt.tsv \
      --chrom "~{chrom}" \
      --out-loci "~{chrom}.trgt_concordance.loci.tsv.gz" \
      --out-summary "~{chrom}.trgt_concordance.summary.json"
  >>>

  output {
    File loci = "~{chrom}.trgt_concordance.loci.tsv.gz"
    File summary = "~{chrom}.trgt_concordance.summary.json"
    File aggregate_summary = "aggregate.summary.json"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
