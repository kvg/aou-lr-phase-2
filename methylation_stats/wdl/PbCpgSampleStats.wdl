version 1.0

# One pb-CpG-tools sample (combined ± hap1/hap2 bedMethyl) → coverage / haplotype stats.
#
# Launch from the Terra aou2_v1_phased_bams data table (one workflow per row):
#   sample_id    = this.aou2_v1_phased_bams_id
#   combined_bed = this.combined_bed
#   hap1_bed     = this.hap1_bed
#   hap2_bed     = this.hap2_bed
#   pbcpg_stats_py = gs://$WORKSPACE_BUCKET/scripts/pbcpg_stats.py
#
# After all rows finish, merge locally or in meth_00_merge_pbcpg_stats.ipynb:
#   python3 scripts/pbcpg_stats.py merge --stats-dir shards --out-dir summaries/manuscript \
#     --covariates covariates.source_rebuilt.csv.gz --discovery-only

workflow PbCpgSampleStats {
  input {
    String sample_id
    File combined_bed
    File? hap1_bed
    File? hap2_bed
    File pbcpg_stats_py

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int cpu = 1
    Int memory_gb = 4
    Int preemptible = 1
    Int disk_gb_floor = 20
    Float disk_gb_multiplier = 2.5
  }

  call SummarizeSample {
    input:
      sample_id = sample_id,
      combined_bed = combined_bed,
      hap1_bed = hap1_bed,
      hap2_bed = hap2_bed,
      pbcpg_stats_py = pbcpg_stats_py,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible = preemptible,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier
  }

  output {
    String sample = SummarizeSample.sample
    File stats_tsv = SummarizeSample.stats_tsv
    File stats_json = SummarizeSample.stats_json
  }

  meta {
    description: "Per-sample pb-CpG-tools bedMethyl coverage and haplotype-track counts."
    allowNestedInputs: true
  }
}

task SummarizeSample {
  input {
    String sample_id
    File combined_bed
    File? hap1_bed
    File? hap2_bed
    File pbcpg_stats_py
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
    Int disk_gb_floor
    Float disk_gb_multiplier
  }

  Float hap_size = size(select_first([hap1_bed, combined_bed]), "GB") + size(select_first([hap2_bed, combined_bed]), "GB")
  Int disk_gb = ceil((size(combined_bed, "GB") + hap_size) * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('sample.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{sample_id}"
    SAMPLE="$(cat sample.txt)"
    echo "[$(date -Is)] pbcpg summarize ${SAMPLE}" >&2

    HAP_ARGS=()
    ~{if defined(hap1_bed) then "HAP_ARGS+=(--hap1 '~{hap1_bed}')" else ""}
    ~{if defined(hap2_bed) then "HAP_ARGS+=(--hap2 '~{hap2_bed}')" else ""}

    python3 "~{pbcpg_stats_py}" summarize \
      --sample-id "${SAMPLE}" \
      --combined "~{combined_bed}" \
      "${HAP_ARGS[@]}" \
      --out-prefix stats

    test -s stats.tsv
    test -s stats.json
  >>>

  output {
    String sample = read_string("sample.txt")
    File stats_tsv = "stats.tsv"
    File stats_json = "stats.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
