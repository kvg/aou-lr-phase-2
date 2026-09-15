version 1.0

# Per-sample Locityper on a streamed CRAM subset.
# Upstream: EichlerLab/AoU_WDL locityper/locityper_stream.wdl (ValidateVariants).
# Prep: notebooks/rw/locityper_00_prep_reference.ipynb
# Submit: notebooks/rw/locityper_01_run_stream.ipynb via `wb workflow`.

workflow LocityperStream {
    input {
        String sample_id

        File cram
        File crai

        File ref_fa_uncompressed
        File ref_fai_uncompressed

        File counts_jf
        File bed

        File locityper_db_tar_gz

        Int locityper_n_cpu = 32
        Int window_grab = 10000
        Int locityper_extra_mem_gb = 6
        Int n_preemptible = 2

        Int max_retry = 5
        Int wait_time = 30

        Int N_split = 200

        # Isaac's Terra run used illumina. PacBio HiFi: "hifi" / "pacbio" (see locityper --help).
        String technology = "illumina"

        # GRCh38 Locityper background interval; appended only for the CRAM subset.
        # Empty string skips the extra grab. CHM13 default is chr17:72950001-77450000.
        String bg_region_bed = "chr17\t72062001\t76562000"

        String locityper_docker = "eichlerlab/locityper:1.4.5.0"
        String util_docker = "python:3.11-slim"
    }

    call SplitBed {
        input:
            bed = bed,
            N = N_split,
            docker = util_docker
    }

    scatter (split_bed in SplitBed.split_beds) {
        call LocityperPreprocessAndGenotype {
            input:
                sample_id = sample_id,
                cram = cram,
                crai = crai,
                reference = ref_fa_uncompressed,
                reference_index = ref_fai_uncompressed,
                db_targz = locityper_db_tar_gz,
                counts_file = counts_jf,
                bed = split_bed,
                locityper_n_cpu = locityper_n_cpu,
                window_grab = window_grab,
                locityper_extra_mem_gb = locityper_extra_mem_gb,
                n_preemptible = n_preemptible,
                max_retry = max_retry,
                wait_time = wait_time,
                technology = technology,
                bg_region_bed = bg_region_bed,
                docker = locityper_docker
        }
    }

    call CombineTarFiles {
        input:
            tar_files = LocityperPreprocessAndGenotype.genotype_tar,
            sample_id = sample_id,
            docker = util_docker
    }

    call Summarize {
        input:
            sample_id = sample_id,
            genotype_tar = CombineTarFiles.combined_tar_gz,
            docker = util_docker
    }

    output {
        File summary_csv = Summarize.summary_csv
        File results_tar_gz = CombineTarFiles.combined_tar_gz
    }
}

task LocityperPreprocessAndGenotype {
    input {
        File cram
        File crai
        File reference
        File reference_index
        File counts_file
        File db_targz
        String sample_id

        File bed

        Int window_grab

        Int locityper_n_cpu

        Int locityper_extra_mem_gb
        Int n_preemptible = 2

        Int max_retry
        Int wait_time

        String technology
        String bg_region_bed
        String docker
    }

    parameter_meta {
        cram: { localization_optional: true }
    }

    Int disk_size = 1 + 4*ceil(size([cram, crai, counts_file, reference, reference_index, db_targz, bed], "GiB"))
    Int locityper_mem_gb = ceil(4.0 * locityper_n_cpu) + locityper_extra_mem_gb

    String output_tar = sample_id + ".locityper.tar.gz"

    command <<<
        set -euo pipefail
        date
        mv ~{reference} reference.fa
        mv ~{reference_index} reference.fa.fai
        mv ~{crai} full.cram.crai

        mkdir -p locityper_preproc

        cp ~{bed} new.bed
        if [ -n "~{bg_region_bed}" ]; then
            echo -e "~{bg_region_bed}" >> new.bed
        fi
        awk -v wgrab="~{window_grab}" '{print $1 "\t" $2 - wgrab "\t" $3 + wgrab "\t" $4}' new.bed > extended.bed
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)

        MAX_RETRIES=~{max_retry}
        CURR_RETRIES=0
        WAIT_TIME=~{wait_time}
        until (( CURR_RETRIES == MAX_RETRIES )) || samtools view -h -T reference.fa  --regions-file extended.bed -C -o subset.cram -X ~{cram} full.cram.crai  >/dev/null 2>&1 ; do
            sleep $WAIT_TIME
            echo $(( CURR_RETRIES++ ))
        done

        if [ ! -f subset.cram ] ; then
            exit 1
        fi

        samtools index subset.cram
        date

        locityper preproc -a subset.cram \
            -r reference.fa \
            -j ~{counts_file} \
            -@ ~{locityper_n_cpu} \
            --technology ~{technology} \
            -o locityper_preproc  >/dev/null 2>&1

        tar -xzf ~{db_targz}

        # Create out_dir and ensure it exists before parallel processing
        mkdir -p out_dir
        mkdir -p out_dir/loci

        process_single_locus() {
            line="$1"
            locus_name=$(echo "$line" | cut -f4)
            # Ensure the locus-specific directory exists
            mkdir -p "out_dir/loci/${locus_name}"

            locityper genotype -a subset.cram \
                -r reference.fa \
                -d vcf_db \
                -p locityper_preproc \
                --subset-loci "${locus_name}" \
                -o out_dir  >/dev/null 2>&1

        }
        export -f process_single_locus

        cat ~{bed} | /usr/bin/parallel --line-buffer -j ~{locityper_n_cpu} process_single_locus {}
        date
        find out_dir -type f -name "*.bam" -exec rm -f {} \;

        rm -f subset.bam
        rm -f new.bed

        tar -czf ~{output_tar} out_dir
        date
    >>>

    output {
        File genotype_tar = output_tar
    }

    runtime {
        memory: "~{locityper_mem_gb} GB"
        cpu: locityper_n_cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: n_preemptible
        docker: docker
    }
}

task SplitBed {
    input {
        File bed
        Int N = 20
        String docker
    }

    Int disk_size = 1 + ceil(size(bed, "GiB"))

    command <<<
        set -euo pipefail

        cat ~{bed} | split -l ~{N} - split_part_ && wc -l split_part_*
    >>>

    output {
        Array[File] split_beds = glob("split_part_*")
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}

task CombineTarFiles {
    input {
        Array[File] tar_files
        String sample_id
        String docker
    }

    Int disk_size = 1 + 10*ceil(size(tar_files, "GiB"))
    String combined_tar = sample_id + ".combined.tar.gz"

    command <<<
        set -euo pipefail

        mkdir -p combined_temp

        for tar_file in ~{sep=" " tar_files}; do
            tar -xzf "$tar_file" -C combined_temp
        done

        tar -czf ~{combined_tar} -C combined_temp .

        rm -rf combined_temp
    >>>

    output {
        File combined_tar_gz = combined_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}

task Summarize {
    input {
        String sample_id
        File genotype_tar
        String docker
    }

    Int disk_size = 1 + 10*ceil(size(genotype_tar, "GiB"))

    command <<<
        set -euo pipefail

        tar -xzvf ~{genotype_tar}

        if [ -d out_dir ]; then
            mv out_dir ~{sample_id}
        elif [ ! -d ~{sample_id} ]; then
            echo "expected out_dir after extracting ~{genotype_tar}" >&2
            exit 1
        fi

        for dir in ~{sample_id}/loci/*/; do
            if [ ! -f "${dir}/res.json.gz" ]; then
                rm -rf "${dir}"
            fi
        done

        python3 - <<'PY'
import gzip
import json
import math
import os
from pathlib import Path

sample_id = "~{sample_id}"
loci_root = Path(sample_id) / "loci"
out_path = Path("gts.filtered.csv")

rows = ["sample\tlocus\tgenotype\tquality\ttotal_reads\tunexpl_reads\tweight_dist\twarnings"]
if loci_root.is_dir():
    for entry in sorted(loci_root.iterdir(), key=lambda p: p.name):
        if not entry.is_dir():
            continue
        json_path = entry / "res.json.gz"
        if not json_path.is_file():
            continue
        with gzip.open(json_path, "rt") as fh:
            res = json.load(fh)
        line = f"{sample_id}\t{entry.name}\t"
        if "genotype" not in res:
            rows.append(line + "*")
            continue
        gt = res["genotype"]
        qual = math.floor(10 * float(res["quality"])) * 0.1
        total_reads = res.get("total_reads", "")
        unexpl_reads = res.get("unexpl_reads", "")
        weight_dist = res.get("weight_dist")
        weight_s = "" if weight_dist is None else f"{float(weight_dist):.5f}"
        warnings = ";".join(res.get("warnings", [])) or "*"
        rows.append(
            line + f"{gt}\t{qual:.1f}\t{total_reads}\t{unexpl_reads}\t{weight_s}\t{warnings}"
        )

out_path.write_text("\n".join(rows) + "\n")
print(f"wrote {len(rows) - 1} loci to {out_path}", flush=True)
PY
    >>>

    output {
        File summary_csv = "gts.filtered.csv"
    }

    runtime {
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}
