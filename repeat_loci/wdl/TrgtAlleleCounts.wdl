version 1.0

# Distinct TRGT allele lengths per TRExplorer locus, for Figure 2C.
# A haplotype's dosage is AL - (END - POS + 1) in bp. See scripts/trgt_allele_counts.py.
#
# Submit twice, once per phase, with the same catalog_bed and a different
# vcf_list and label. vcf_list is one TRGT VCF URI per line, or
# sample_id<TAB>URI with no header (the same shape as trgt_table.txt).
# Phase 2 can reuse the TrgtPlvi list. Phase 1 is trgt.phase1.txt.

workflow TrgtAlleleCounts {
  input {
    File vcf_list
    File catalog_bed
    Int shard_size = 50
    String label = "trgt_alleles"

    File trgt_allele_counts_py

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int cpu = 4
    Int memory_gb = 16
    Int preemptible = 2
  }

  call SplitList {
    input:
      vcf_list = vcf_list,
      shard_size = shard_size,
      docker = docker
  }

  scatter (shard in SplitList.shards) {
    Array[File] shard_vcfs = read_lines(shard)
    call ShardCounts {
      input:
        vcfs = shard_vcfs,
        catalog_bed = catalog_bed,
        trgt_allele_counts_py = trgt_allele_counts_py,
        docker = docker,
        cpu = cpu,
        memory_gb = memory_gb,
        preemptible = preemptible
    }
  }

  call MergeAlleles {
    input:
      counts = ShardCounts.counts,
      catalog_bed = catalog_bed,
      label = label,
      trgt_allele_counts_py = trgt_allele_counts_py,
      docker = docker
  }

  output {
    File alleles = MergeAlleles.alleles
    File summary = MergeAlleles.summary
    Array[File] shard_counts = ShardCounts.counts
  }

  meta {
    description: "TRGT allele-length counts per catalog locus (Figure 2C rarefaction and example histograms)."
    allowNestedInputs: true
  }
}

task SplitList {
  input {
    File vcf_list
    Int shard_size
    String docker
  }

  command <<<
    set -euo pipefail
    mkdir shards
    # URI-only lines pass through. sample_id<TAB>URI lines (Phase 1, trgt_table.txt) keep the path.
    grep -v '^[[:space:]]*$' "~{vcf_list}" | awk -F'\t' 'NF > 1 { print $2; next } { print $1 }' | split -l ~{shard_size} -d -a 4 - shards/shard_
    ls shards | wc -l
  >>>

  output {
    Array[File] shards = glob("shards/shard_*")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}

task ShardCounts {
  input {
    Array[File] vcfs
    File catalog_bed
    File trgt_allele_counts_py
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
  }

  Int disk_gb = ceil(size(vcfs, "GB") * 1.2) + 20

  command <<<
    set -euo pipefail
    python3 "~{trgt_allele_counts_py}" counts \
      --catalog-bed "~{catalog_bed}" \
      --vcf-list "~{write_lines(vcfs)}" \
      --threads ~{cpu} \
      --out counts.npz
  >>>

  output {
    File counts = "counts.npz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task MergeAlleles {
  input {
    Array[File] counts
    File catalog_bed
    String label
    File trgt_allele_counts_py
    String docker
  }

  Int disk_gb = ceil(size(counts, "GB") * 2) + 20

  command <<<
    set -euo pipefail
    python3 "~{trgt_allele_counts_py}" merge \
      --catalog-bed "~{catalog_bed}" \
      --counts ~{sep=" " counts} \
      --out "~{label}.alleles.tsv.gz" \
      --out-summary "~{label}.summary.json"
  >>>

  output {
    File alleles = "~{label}.alleles.tsv.gz"
    File summary = "~{label}.summary.json"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
  }
}
