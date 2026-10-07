version 1.0

# Per-locus PLVI (SD of longest-pure-segment length across haplotypes;
# Danzi et al. 2025) from per-sample TRGT VCFs, ranked within motif-length
# groups. Used to choose the Figure 2C example loci. See scripts/trgt_plvi.py.
#
# vcf_list: one TRGT VCF URI per line (a random participant subsample is
# enough; LPS stats are summed across shards).

workflow TrgtPlvi {
  input {
    File vcf_list
    File catalog_bed
    Int shard_size = 50
    Float min_call_rate = 0.8
    String label = "trgt_plvi"

    File trgt_plvi_py

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
    call ShardStats {
      input:
        vcfs = shard_vcfs,
        catalog_bed = catalog_bed,
        trgt_plvi_py = trgt_plvi_py,
        docker = docker,
        cpu = cpu,
        memory_gb = memory_gb,
        preemptible = preemptible
    }
  }

  call MergePlvi {
    input:
      stats = ShardStats.stats,
      min_call_rate = min_call_rate,
      label = label,
      trgt_plvi_py = trgt_plvi_py,
      docker = docker
  }

  output {
    File plvi = MergePlvi.plvi
    Array[File] shard_stats = ShardStats.stats
  }

  meta {
    description: "PLVI per TRGT locus from a TRGT VCF subsample, ranked within motif-length groups."
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
    grep -v '^[[:space:]]*$' "~{vcf_list}" | split -l ~{shard_size} -d -a 4 - shards/shard_
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

task ShardStats {
  input {
    Array[File] vcfs
    File catalog_bed
    File trgt_plvi_py
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
  }

  Int disk_gb = ceil(size(vcfs, "GB") * 1.2) + 20

  command <<<
    set -euo pipefail
    python3 "~{trgt_plvi_py}" stats \
      --catalog-bed "~{catalog_bed}" \
      --vcf-list "~{write_lines(vcfs)}" \
      --threads ~{cpu} \
      --out stats.tsv.gz
  >>>

  output {
    File stats = "stats.tsv.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task MergePlvi {
  input {
    Array[File] stats
    Float min_call_rate
    String label
    File trgt_plvi_py
    String docker
  }

  Int disk_gb = ceil(size(stats, "GB") * 2) + 20

  command <<<
    set -euo pipefail
    python3 "~{trgt_plvi_py}" merge \
      --stats ~{sep=" " stats} \
      --min-call-rate ~{min_call_rate} \
      --out "~{label}.tsv.gz"
  >>>

  output {
    File plvi = "~{label}.tsv.gz"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
  }
}
