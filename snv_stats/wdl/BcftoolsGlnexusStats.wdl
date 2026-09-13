version 1.0

# One DeepVariant + GLnexus chromosome shard → bcftools stats.
#
# Launch from the Terra GL_INTERVAL_set data table (one workflow per row):
#   chrom = this.GL_INTERVAL_set_id
#   vcf   = this.VCF
#
# After chr1–22 / X / Y finish, merge locally:
#   python3 scripts/snv_bcftools_sample_qc.py --merge-dir shards --out-dir summaries

workflow BcftoolsGlnexusStats {
  input {
    String chrom
    File vcf

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    # GLnexus typically leaves FILTER as "." (unfiltered), not PASS.
    # bcftools -f PASS then yields 0 records. Include both.
    String apply_filters = "PASS,."
    Int cpu = 4
    Int memory_gb = 8
    Int preemptible = 1
    Int disk_gb_floor = 80
    Float disk_gb_multiplier = 2.0
    Int bcftools_threads = 4
  }

  call StatsShard {
    input:
      chrom = chrom,
      vcf = vcf,
      apply_filters = apply_filters,
      bcftools_threads = bcftools_threads,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier,
      preemptible = preemptible
  }

  output {
    String chrom_id = StatsShard.chrom_id
    File stats = StatsShard.stats
  }

  meta {
    description: "bcftools stats -s - on one GLnexus chrom shard (Terra data-table row)."
    allowNestedInputs: true
  }
}

task StatsShard {
  input {
    String chrom
    File vcf
    String apply_filters
    Int bcftools_threads
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb_floor
    Float disk_gb_multiplier
    Int preemptible
  }

  Int disk_gb = ceil(size(vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('chrom.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{chrom}"
    CHROM="$(cat chrom.txt)"
    echo "[$(date -Is)] bcftools stats ${CHROM}" >&2
    ls -lh "~{vcf}" || true
    df -h . || true

    FILTER_ARGS=()
    if [[ -n "~{apply_filters}" ]]; then
      FILTER_ARGS=(-f "~{apply_filters}")
    fi

    # Cheap peek: FILTER on the first sites *without* -f, so an empty stats
    # file is diagnosable (GLnexus "." vs PASS).
    set +o pipefail
    bcftools view -H -G --threads 1 "~{vcf}" \
      | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' \
      | head -n 20 > sites_head.tsv || true
    set -o pipefail
    echo "=== first sites (CHROM POS REF ALT FILTER) ===" >&2
    cat sites_head.tsv >&2 || true

    bcftools stats -s - --threads ~{bcftools_threads} "${FILTER_ARGS[@]}" "~{vcf}" \
      > "${CHROM}.stats.txt"
    ln -s "${CHROM}.stats.txt" stats.txt

    python3 - <<'PY'
from pathlib import Path
p = Path("stats.txt")
text = p.read_text(errors="replace")
if "# PSC" not in text or "\nPSC\t" not in text or not text.endswith("\n"):
    raise SystemExit(f"incomplete bcftools stats: {p} ({p.stat().st_size} bytes)")
n_records = 0
for line in text.splitlines():
    if line.startswith("SN\t") and "number of records:" in line:
        n_records = int(line.rsplit("\t", 1)[-1])
        break
print(f"OK {p.resolve()}: {p.stat().st_size} bytes, records={n_records}")
if n_records == 0:
    head = Path("sites_head.tsv").read_text() if Path("sites_head.tsv").exists() else ""
    print(
        "WARNING: 0 records after apply-filters. If sites_head shows FILTER=., "
        "use apply_filters=PASS,. (GLnexus default is unfiltered '.' not PASS).",
        flush=True,
    )
    print(head, flush=True)
PY
  >>>

  output {
    String chrom_id = read_string("chrom.txt")
    File stats = "stats.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
