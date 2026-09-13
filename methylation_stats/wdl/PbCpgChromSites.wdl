version 1.0

# Extract one chromosome of a combined bedMethyl (coverage-filtered) for
# Primrose vs Jasmine concordance. Do not launch on all 12k rows.
# Import summaries/manuscript/pbcpg_concordance.entities.tsv as table
# `pbcpg_concordance` (25 Primrose + 25 Jasmine), then:
#
#   sample_id    = this.pbcpg_concordance_id
#   combined_bed = this.combined_bed
#   pbcpg_stats_py = gs://$WORKSPACE_BUCKET/scripts/pbcpg_stats.py
#
# Write output `sites` back onto that table. meth_00 then pull-dumps.

workflow PbCpgChromSites {
  input {
    String sample_id
    File combined_bed
    File pbcpg_stats_py
    String chrom = "chr22"
    Int min_cov = 10

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int cpu = 1
    Int memory_gb = 4
    Int preemptible = 1
    Int disk_gb_floor = 20
    Float disk_gb_multiplier = 2.5
  }

  call ExtractChrom {
    input:
      sample_id = sample_id,
      combined_bed = combined_bed,
      pbcpg_stats_py = pbcpg_stats_py,
      chrom = chrom,
      min_cov = min_cov,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible = preemptible,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier
  }

  output {
    String sample = ExtractChrom.sample
    File sites = ExtractChrom.sites
  }

  meta {
    description: "Chromosome-restricted pb-CpG-tools site dump for caller concordance."
    allowNestedInputs: true
  }
}

task ExtractChrom {
  input {
    String sample_id
    File combined_bed
    File pbcpg_stats_py
    String chrom
    Int min_cov
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
    Int disk_gb_floor
    Float disk_gb_multiplier
  }

  Int disk_gb = ceil(size(combined_bed, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('sample.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{sample_id}"
    SAMPLE="$(cat sample.txt)"
    python3 -c "import pathlib,sys; pathlib.Path('chrom.txt').write_text(sys.argv[1].strip() + chr(10))" "~{chrom}"
    CHROM="$(cat chrom.txt)"

    python3 "~{pbcpg_stats_py}" extract-chrom \
      --bed "~{combined_bed}" \
      --chrom "${CHROM}" \
      --min-cov ~{min_cov} \
      --sample-id "${SAMPLE}" \
      --out "${SAMPLE}.${CHROM}.sites.tsv.gz"

    ln -s "${SAMPLE}.${CHROM}.sites.tsv.gz" sites.tsv.gz
  >>>

  output {
    String sample = read_string("sample.txt")
    File sites = "sites.tsv.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
