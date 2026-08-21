version 1.0

# Genome-wide Tractor-Mix for Terra Methods Repository (self-contained; no imports).
#
# Shape:
#   MakeGRM once (separate grm_vcfs; keep small — do not dump all autosomes here)
#   FitNull once per phenotype
#   Extract once per chromosome
#   Score once per (phenotype × chromosome)
#
# Terra chrom-set launch: pass parallel Arrays
#   flare_vcfs <- this.chrom_set.model_chr_anc_vcf
#   chroms     <- this.chrom_set.aou_lr_chrom_id  (or explicit chr1..chr22 names)
# in matching order. generate URIs with scripts/resolve_flare_uris.py --autosomes.

task CheckChromVcfPairs {
  input {
    # Lengths only — do NOT pass Array[File] here or Cromwell localizes every
    # FLARE VCF onto this tiny disk (chr1 alone OOMs a 10 GB boot disk).
    Int n_chroms
    Int n_vcfs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
  }

  command <<<
    set -euo pipefail
    N_CHROM=~{n_chroms}
    N_VCF=~{n_vcfs}
    if [[ "${N_CHROM}" -ne "${N_VCF}" ]]; then
      echo "chroms length (${N_CHROM}) != flare_vcfs length (${N_VCF})" >&2
      exit 1
    fi
    if [[ "${N_CHROM}" -lt 1 ]]; then
      echo "Need at least one chromosome / FLARE VCF" >&2
      exit 1
    fi
    echo "OK: ${N_CHROM} chrom × VCF pairs"
  >>>

  output {
    Int n_chroms_out = n_chroms
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task ExtractTractsFlare {
  input {
    File flare_vcf
    String chrom
    File analysis_samples
    Int num_ancs
    File reorder_dosages_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
    Int cpu = 4
    Int memory_gb = 16
    # Floor added on top of size-based disk (see runtime). Override via workflow.
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.5
    Int preemptible = 0
  }

  # Localized VCF + extract outputs (dosages) + headroom. chr1 FLARE is huge.
  Int disk_gb = ceil(size(flare_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p extract ordered

    echo "chrom=~{chrom}"
    echo "Host memory (MB):"; free -m || true
    echo "Disk:"; df -h . || true
    echo "Starting extract-tracts-flare at $(date -Is)"

    extract-tracts-flare \
      --vcf "~{flare_vcf}" \
      --num-ancs ~{num_ancs} \
      --output-dir extract \
      --compress-output \
      --samples "~{analysis_samples}"

    echo "Finished extract-tracts-flare at $(date -Is)"
    ls -lh extract/ || true

    python3 - <<'PY'
import re
import shutil
from pathlib import Path

src = Path("extract")
od = Path("ordered")
for p in src.glob("*.dosage.txt.gz"):
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if not m:
        raise SystemExit(f"unexpected dosage name: {p.name}")
    shutil.copy2(p, od / f"anc_{int(m.group(1)):02d}.dosage.txt.gz")
for p in src.glob("*.hapcount.txt.gz"):
    m = re.search(r"\.anc(\d+)\.hapcount\.txt\.gz$", p.name)
    if not m:
        raise SystemExit(f"unexpected hapcount name: {p.name}")
    shutil.copy2(p, od / f"anc_{int(m.group(1)):02d}.hapcount.txt.gz")
order = src / "dosage_sample_order.txt"
if not order.exists():
    raise SystemExit("missing extract/dosage_sample_order.txt")
shutil.copy2(order, od / "dosage_sample_order.txt")
PY
    cp "~{analysis_samples}" ordered/analysis_samples.txt
    # Strip accidental JSON quotes from Terra Array[String] members.
    python3 -c "import pathlib,sys; pathlib.Path('ordered/chrom.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{chrom}"

    python3 - <<'PY'
from pathlib import Path
import gzip

chrom = Path("ordered/chrom.txt").read_text().strip()
samples = Path("ordered/analysis_samples.txt").read_text().split()
dfiles = sorted(Path("ordered").glob("anc_*.dosage.txt.gz"))
assert dfiles, "no ordered dosage files"
with gzip.open(dfiles[0], "rt") as fh:
    got = fh.readline().rstrip("\n").split("\t")[5:]
if got != samples:
    raise SystemExit(
        f"sample order mismatch: dosage has {len(got)}, analysis has {len(samples)}"
    )
print(f"OK {chrom}: {len(samples)} samples, {len(dfiles)} dosage files")
PY
  >>>

  output {
    String chrom_id = read_string("ordered/chrom.txt")
    Array[File] dosage_files = glob("ordered/anc_*.dosage.txt.gz")
    Array[File] hapcount_files = glob("ordered/anc_*.hapcount.txt.gz")
    File dosage_sample_order = "ordered/dosage_sample_order.txt"
    File chrom_txt = "ordered/chrom.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task MakeSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File sparsify_grm_script
    File make_plink_keep_script
    Float kinship_threshold = 0.05
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb_floor = 100
    # reheader writes a full copy; concat may write another; bed/rel on top.
    Float disk_gb_multiplier = 4.0
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(grm_vcfs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p grm

    N_VCF=0
    for vcf in ~{sep=" " grm_vcfs}; do
      echo "$vcf" >> grm/vcf_list.txt
      N_VCF=$((N_VCF+1))
    done

    if [[ "$N_VCF" -eq 1 ]]; then
      INPUT_VCF=$(head -n1 grm/vcf_list.txt)
    else
      bcftools concat -f grm/vcf_list.txt -Oz -o grm/concat.vcf.gz --threads ~{cpu}
      bcftools index -t grm/concat.vcf.gz
      INPUT_VCF=grm/concat.vcf.gz
    fi

    # FLARE anc VCFs often declare ##FORMAT=<ID=GT,...> twice. Recent PLINK2
    # errors with "Duplicate FORMAT/GT header line" (pilot chr1+chr22 concat
    # sometimes masked this; single-VCF GRM hits the raw header).
    bcftools view -h "$INPUT_VCF" \
      | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
      > grm/header.fixed.txt
    echo "FORMAT/GT header lines: $(bcftools view -h "$INPUT_VCF" | grep -c '^##FORMAT=<ID=GT,' || true) -> $(grep -c '^##FORMAT=<ID=GT,' grm/header.fixed.txt || true)"
    bcftools reheader -h grm/header.fixed.txt -o grm/plink_in.vcf.gz "$INPUT_VCF"
    bcftools index -t grm/plink_in.vcf.gz || bcftools index -c grm/plink_in.vcf.gz
    INPUT_VCF=grm/plink_in.vcf.gz

    echo "=== analysis_samples diagnostics ==="
    ls -l "~{analysis_samples}" || true
    wc -l "~{analysis_samples}" || true
    head -n 5 "~{analysis_samples}" || true
    echo "=== end diagnostics ==="

    bcftools query -l "$INPUT_VCF" > grm/vcf_samples.txt
    python3 "~{make_plink_keep_script}" \
      --analysis-samples "~{analysis_samples}" \
      --fam grm/vcf_samples.txt \
      --out-keep grm/keep_precheck.txt \
      --out-samples grm/analysis_samples.intersect.txt
    echo "=== keep_precheck head ==="
    head -n 5 grm/keep_precheck.txt || true

    plink2 --vcf "$INPUT_VCF" \
      --double-id \
      --keep grm/keep_precheck.txt \
      --make-bed \
      --set-all-var-ids "@:#:\$r:\$a" \
      --out grm/all \
      --threads ~{cpu}

    python3 "~{make_plink_keep_script}" \
      --analysis-samples grm/analysis_samples.intersect.txt \
      --fam grm/all.fam \
      --out-keep grm/keep_iid.txt \
      --out-samples grm/analysis_samples.intersect.txt
    echo "=== keep_iid / fam head ==="
    head -n 3 grm/keep_iid.txt grm/all.fam || true

    plink2 --bfile grm/all \
      --keep grm/keep_iid.txt \
      --indep-pairwise 500kb 0.2 \
      --out grm/prune \
      --threads ~{cpu}

    plink2 --bfile grm/all \
      --keep grm/keep_iid.txt \
      --extract grm/prune.prune.in \
      --make-rel square \
      --out grm/pruned \
      --threads ~{cpu}

    Rscript "~{sparsify_grm_script}" \
      --rel grm/pruned.rel \
      --rel-id grm/pruned.rel.id \
      --samples grm/analysis_samples.intersect.txt \
      --threshold ~{kinship_threshold} \
      --out-rds grm_sparse.rds \
      --out-counts grm_relatedness_summary.tsv \
      --out-histogram grm_kinship_histogram.tsv \
      --out-bands grm_relationship_bands.tsv \
      --out-close-pairs grm_close_pairs.tsv

    cp grm/pruned.rel .
    cp grm/pruned.rel.id .
    cp grm/keep_iid.txt plink_keep_iid.txt
    cp grm/analysis_samples.intersect.txt analysis_samples.intersect.txt
  >>>

  output {
    File grm_sparse_rds = "grm_sparse.rds"
    File grm_summary = "grm_relatedness_summary.tsv"
    File grm_histogram = "grm_kinship_histogram.tsv"
    File grm_bands = "grm_relationship_bands.tsv"
    File grm_close_pairs = "grm_close_pairs.tsv"
    File rel_matrix = "pruned.rel"
    File rel_ids = "pruned.rel.id"
    File plink_keep = "plink_keep_iid.txt"
    File analysis_samples_intersect = "analysis_samples.intersect.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task FitNull {
  input {
    File pheno_cov
    String phenotype
    File covariate_columns
    File grm_sparse_rds
    File fit_null_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 50
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail

    Rscript "~{fit_null_script}" \
      --pheno-cov "~{pheno_cov}" \
      --phenotype "~{phenotype}" \
      --covariates "~{covariate_columns}" \
      --grm-rds "~{grm_sparse_rds}" \
      --out-null-rds "~{phenotype}.null_model.rds" \
      --out-null-export null_export

    tar czf null_export.tar.gz null_export
  >>>

  output {
    File null_model_rds = "~{phenotype}.null_model.rds"
    File null_export_tar = "null_export.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task Score {
  input {
    String phenotype
    String chrom
    File null_export_tar
    Array[File] dosage_files
    Int ac_threshold = 50
    Int score_threads = 8
    Int chunk_size = 2048
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 200
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    # Sanitize in case Terra Array[String] members include JSON quotes.
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    OUT="${PHENO}.${CHROM}.tractor_mix.tsv"

    tar xzf "~{null_export_tar}"
    DOSAGE_ARGS=()
    for f in ~{sep=" " dosage_files}; do
      DOSAGE_ARGS+=("$f")
    done

    tractor-mix-score \
      --null-export null_export \
      --dosage-files "${DOSAGE_ARGS[@]}" \
      --out "${OUT}" \
      --ac-threshold ~{ac_threshold} \
      --threads ~{score_threads} \
      --chunk-size ~{chunk_size}

    # Stable name for WDL output declaration (avoids quote chars in ~{chrom}).
    cp "${OUT}" results.tractor_mix.tsv
    printf '%s\n' "${OUT}" > results.output_name.txt
  >>>

  output {
    File results_tsv = "results.tractor_mix.tsv"
    File results_name = "results.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

workflow TractorMixGenome {
  input {
    # Parallel arrays (same length / order). Terra chrom-set:
    #   flare_vcfs = this.SET.model_chr_anc_vcf
    #   chroms     = explicit chr names matching those files
    Array[File] flare_vcfs
    Array[String] chroms
    Int num_ancs = 5

    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns

    # Keep GRM chroms separate and usually small (e.g. chr1+chr22).
    Array[File] grm_vcfs

    File reorder_dosages_script
    File sparsify_grm_script
    File fit_null_script
    File make_plink_keep_script

    Float kinship_threshold = 0.05
    Int ac_threshold = 50
    Int score_threads = 8
    Int chunk_size = 2048

    Int extract_cpu = 4
    Int extract_memory_gb = 16
    # Added to size(flare_vcf)*multiplier for Extract (not a fixed total disk).
    Int extract_disk_gb_floor = 100
    Float extract_disk_gb_multiplier = 3.5
    Int score_disk_gb = 200
    Int grm_disk_gb_floor = 100
    Float grm_disk_gb_multiplier = 4.0

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)

  call CheckChromVcfPairs as Check {
    input:
      n_chroms = length(chroms),
      n_vcfs = length(flare_vcfs),
      docker = docker
  }

  call MakeSparseGRM as MakeGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      sparsify_grm_script = sparsify_grm_script,
      make_plink_keep_script = make_plink_keep_script,
      kinship_threshold = kinship_threshold,
      docker = docker,
      disk_gb_floor = grm_disk_gb_floor,
      disk_gb_multiplier = grm_disk_gb_multiplier
  }

  scatter (pheno in phenotypes) {
    call FitNull {
      input:
        pheno_cov = pheno_cov,
        phenotype = pheno,
        covariate_columns = covariate_columns,
        grm_sparse_rds = MakeGRM.grm_sparse_rds,
        fit_null_script = fit_null_script,
        docker = docker
    }
  }

  scatter (pair in zip(chroms, flare_vcfs)) {
    String chrom = pair.left
    File flare_vcf = pair.right

    call ExtractTractsFlare as Extract {
      input:
        flare_vcf = flare_vcf,
        chrom = chrom,
        analysis_samples = analysis_samples,
        num_ancs = num_ancs,
        reorder_dosages_script = reorder_dosages_script,
        docker = docker,
        cpu = extract_cpu,
        memory_gb = extract_memory_gb,
        disk_gb_floor = extract_disk_gb_floor,
        disk_gb_multiplier = extract_disk_gb_multiplier
    }

    scatter (i in range(length(phenotypes))) {
      call Score {
        input:
          phenotype = phenotypes[i],
          chrom = chrom,
          null_export_tar = FitNull.null_export_tar[i],
          dosage_files = Extract.dosage_files,
          ac_threshold = ac_threshold,
          score_threads = score_threads,
          chunk_size = chunk_size,
          docker = docker,
          disk_gb = score_disk_gb
      }
    }
  }

  output {
    Int n_chroms = Check.n_chroms_out
    File grm_sparse_rds = MakeGRM.grm_sparse_rds
    File grm_summary = MakeGRM.grm_summary
    File grm_histogram = MakeGRM.grm_histogram
    File grm_bands = MakeGRM.grm_bands
    File grm_close_pairs = MakeGRM.grm_close_pairs
    Array[File] null_model_rds = FitNull.null_model_rds
    # Outer = chrom, inner = phenotype (same order as chroms / phenotypes).
    Array[Array[File]] results_tsvs = Score.results_tsv
    Array[Array[File]] dosage_files = Extract.dosage_files
    Array[String] chrom_ids = Extract.chrom_id
  }

  meta {
    description: "Genome-wide Tractor-Mix: shared sparse GRM + nulls; per-chr FLARE extract and score."
    allowNestedInputs: true
  }
}
