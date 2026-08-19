version 1.0

# This file is intentionally self-contained so it can be imported directly
# into a Terra Methods Repository, which does not resolve sibling WDL imports.

task ExtractTractsFlare {
  input {
    File flare_vcf
    File analysis_samples
    Int num_ancs
    # Kept for Terra config compatibility; sample reorder is done in-extract.
    File reorder_dosages_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb = 500
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p extract ordered

    echo "Host memory (MB):"; free -m || true
    echo "Disk:"; df -h . || true
    echo "Starting extract-tracts-flare at $(date -Is)"
    ls -l /usr/local/bin/extract-tracts-flare
    extract-tracts-flare --help | head -n 5

    extract-tracts-flare \
      --vcf "~{flare_vcf}" \
      --num-ancs ~{num_ancs} \
      --output-dir extract \
      --compress-output \
      --samples "~{analysis_samples}"

    echo "Finished extract-tracts-flare at $(date -Is)"
    ls -lh extract/ || true
    df -h . || true

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

    python3 - <<'PY'
from pathlib import Path
import gzip

samples = Path("~{analysis_samples}").read_text().split()
dfiles = sorted(Path("ordered").glob("anc_*.dosage.txt.gz"))
assert dfiles, "no ordered dosage files"
with gzip.open(dfiles[0], "rt") as fh:
    got = fh.readline().rstrip("\n").split("\t")[5:]
if got != samples:
    raise SystemExit(f"sample order mismatch: dosage has {len(got)}, analysis has {len(samples)}")
print(f"OK: {len(samples)} samples ordered across {len(dfiles)} dosage files")
PY
  >>>

  output {
    Array[File] dosage_files = glob("ordered/anc_*.dosage.txt.gz")
    Array[File] hapcount_files = glob("ordered/anc_*.hapcount.txt.gz")
    File dosage_sample_order = "ordered/dosage_sample_order.txt"
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
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 2
  }

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

    # Quick ID overlap check against VCF before the expensive PLINK convert
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

    # Rebuild keep against the FAM (FID=IID pairs for --double-id beds)
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

task FitNullAndScore {
  input {
    File pheno_cov
    String phenotype
    File covariate_columns
    File grm_sparse_rds
    Array[File] dosage_files
    File fit_null_and_score_script
    Int ac_threshold = 50
    Int n_core = 8
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
    Int cpu = 8
    # 16 GB OOM-killed TractorMix.score after glmmkin (5 ancestry dosages ×
    # ~9k samples, foreach workers). 64 GB leaves headroom if Sigma_i is dense.
    Int memory_gb = 64
    Int disk_gb = 100
    # Score is a multi-hour uncheckpointed R process; preemptible VMs waste
    # the whole null-model + TractorMix.score run when reclaimed.
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    DOSAGE_ARGS=()
    for f in ~{sep=" " dosage_files}; do
      DOSAGE_ARGS+=("$f")
    done

    Rscript "~{fit_null_and_score_script}" \
      --pheno-cov "~{pheno_cov}" \
      --phenotype "~{phenotype}" \
      --covariates "~{covariate_columns}" \
      --grm-rds "~{grm_sparse_rds}" \
      --dosage-files "${DOSAGE_ARGS[@]}" \
      --tractor-mix-score-r "/opt/Tractor-Mix/TractorMix.score.R" \
      --ac-threshold ~{ac_threshold} \
      --n-core ~{n_core} \
      --out-tsv "~{phenotype}.tractor_mix.tsv" \
      --out-null-rds "~{phenotype}.null_model.rds"
  >>>

  output {
    File results_tsv = "~{phenotype}.tractor_mix.tsv"
    File null_model_rds = "~{phenotype}.null_model.rds"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

workflow TractorMixPilot {
  input {
    # chr22 (or other) phased VCF with FLARE AN1/AN2 annotations
    File flare_vcf
    Int num_ancs

    # Sample / phenotype table from notebooks/tractor_01_prepare_inputs.ipynb.
    # Use the shared recommended-full-complete cohort for all matched runs.
    # `pheno_cov` contains phenotypes and both covariate matrices; select one
    # with covariate_columns_limited.txt or covariate_columns_full.txt.
    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns

    # VCFs used to build the GRM (chr22-only OK for pilot; more chroms better)
    Array[File] grm_vcfs

    # Pipeline scripts (from this repo)
    File reorder_dosages_script
    File sparsify_grm_script
    File fit_null_and_score_script
    File make_plink_keep_script

    Float kinship_threshold = 0.05
    Int ac_threshold = 50
    Int score_n_core = 8

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.3.2"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)

  call ExtractTractsFlare as Extract {
    input:
      flare_vcf = flare_vcf,
      analysis_samples = analysis_samples,
      num_ancs = num_ancs,
      reorder_dosages_script = reorder_dosages_script,
      docker = docker
  }

  call MakeSparseGRM as MakeGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      sparsify_grm_script = sparsify_grm_script,
      make_plink_keep_script = make_plink_keep_script,
      kinship_threshold = kinship_threshold,
      docker = docker
  }

  scatter (pheno in phenotypes) {
    call FitNullAndScore as Score {
      input:
        pheno_cov = pheno_cov,
        phenotype = pheno,
        covariate_columns = covariate_columns,
        grm_sparse_rds = MakeGRM.grm_sparse_rds,
        dosage_files = Extract.dosage_files,
        fit_null_and_score_script = fit_null_and_score_script,
        ac_threshold = ac_threshold,
        n_core = score_n_core,
        docker = docker
    }
  }

  output {
    Array[File] dosage_files = Extract.dosage_files
    Array[File] hapcount_files = Extract.hapcount_files
    File grm_sparse_rds = MakeGRM.grm_sparse_rds
    File grm_summary = MakeGRM.grm_summary
    File grm_histogram = MakeGRM.grm_histogram
    File grm_bands = MakeGRM.grm_bands
    File grm_close_pairs = MakeGRM.grm_close_pairs
    Array[File] results_tsvs = Score.results_tsv
    Array[File] null_model_rds = Score.null_model_rds
  }

  meta {
    description: "Tractor-Mix chr22 pilot: Rust FLARE extract tracts, sparse GRM, per-phenotype score tests."
    allowNestedInputs: true
  }
}
