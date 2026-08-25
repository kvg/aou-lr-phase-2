version 1.0

# Per-chromosome union of SNV/indel + SV + FLARE into one annotated BCF.
#
# Pipeline:
#   1) PrepareChrom  – chrom-subset SV if needed; index SNV/FLARE (no full sort)
#   2) MakeRegions   – split chromosome into ~equal-variant shards (SNV density)
#   3) scatter PropagateShard – region-subset + streaming merge → shard BCF
#   4) ConcatShards  – bcftools concat + merged sanity-check stats
#
# Sparse FORMAT (default): each site only emits FORMAT tags from matched sources.
# Rayon --threads parallelizes sample assembly within each shard.

task PrepareChrom {
  input {
    String chromosome
    File? snv_indel_vcf
    File? sv_vcf
    File? sv_vcf_index  # optional .csi; enables fast -r extract
    File? flare_vcf
    Boolean sv_is_whole_genome = true
    String docker
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 2.0
    Int preemptible = 0
  }

  Float snv_gb = if defined(snv_indel_vcf) then size(select_first([snv_indel_vcf]), "GB") else 0.0
  Float sv_gb = if defined(sv_vcf) then size(select_first([sv_vcf]), "GB") else 0.0
  Float flare_gb = if defined(flare_vcf) then size(select_first([flare_vcf]), "GB") else 0.0
  # Localized inputs + chrom extracts (no full-file sort temp).
  Int disk_gb = ceil((snv_gb + sv_gb + flare_gb) * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM="~{chromosome}"
    echo "[$(date -Is)] PrepareChrom ${CHROM}" >&2
    df -h . || true

    # Place an index beside a VCF/BCF so htslib can find it (Cromwell often
    # localizes File? indexes to a path that is not adjacent to the callset).
    attach_index() {
      local src="$1"
      local idx="$2"
      if [ -z "${idx}" ] || [ ! -f "${idx}" ]; then
        return 1
      fi
      # Prefer the conventional "<file>.csi" name next to src.
      ln -sfn "$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${idx}")" "${src}.csi"
      return 0
    }

    prepare_vcf() {
      local src="$1"
      local dst="$2"
      local filter_chrom="$3"
      local optional_index="${4:-}"

      # Per-chromosome SNV/FLARE are already sorted callsets. Copy + index only.
      # A full bcftools sort on ~3–4M sites × ~12k samples is multi-hour and
      # was the dominant cost in PrepareChrom.
      if [ "${filter_chrom}" != "true" ]; then
        echo "Link/index (no sort): $(basename "${src}") -> ${dst}" >&2
        if [[ "${src}" == *.vcf.gz || "${src}" == *.vcf.bgz ]]; then
          cp "${src}" "${dst}"
        elif [[ "${src}" == *.bcf || "${src}" == *.bcf.gz ]]; then
          bcftools view -Oz --threads ~{cpu} -o "${dst}" "${src}"
        else
          bcftools view -Oz --threads ~{cpu} -o "${dst}" "${src}"
        fi
        rm -f "${dst}.tbi" "${dst}.csi"
        if ! bcftools index -t "${dst}" 2>/dev/null; then
          bcftools index -c "${dst}"
        fi
        echo -n "$(basename "${dst}") records: "
        bcftools index -n "${dst}" 2>/dev/null || bcftools view -H "${dst}" | wc -l
        return 0
      fi

      # Whole-genome → chromosome extract. Prefer indexed -r; fall back to -t.
      # Extraction preserves POS order; same-POS REF/ALT order is normalized in
      # the Rust merger (bcftools does not sort by REF/ALT within a position).
      local view_region_args=()
      if attach_index "${src}" "${optional_index}" \
        || [ -f "${src}.csi" ] || [ -f "${src}.tbi" ] || [ -f "${src}.crai" ]; then
        view_region_args=(-r "${CHROM}")
        echo "Extract ${CHROM} via -r (indexed): $(basename "${src}") -> ${dst}" >&2
      else
        view_region_args=(-t "${CHROM}")
        echo "Extract ${CHROM} via -t (no index; streaming whole file): $(basename "${src}") -> ${dst}" >&2
      fi

      bcftools view "${view_region_args[@]}" -Oz --threads ~{cpu} -o "${dst}" "${src}"
      rm -f "${dst}.tbi" "${dst}.csi"
      if ! bcftools index -t "${dst}" 2>/dev/null; then
        bcftools index -c "${dst}"
      fi
      echo -n "$(basename "${dst}") records: "
      bcftools index -n "${dst}" 2>/dev/null || bcftools view -H "${dst}" | wc -l
    }

    HAS_SNV=0
    HAS_SV=0
    HAS_FLARE=0

    if [ -n "~{snv_indel_vcf}" ] && [ "~{snv_indel_vcf}" != "" ]; then
      prepare_vcf "~{snv_indel_vcf}" snv_indel.chr.vcf.gz false
      HAS_SNV=1
    fi
    if [ -n "~{sv_vcf}" ] && [ "~{sv_vcf}" != "" ]; then
      SV_IDX=""
      if [ -n "~{sv_vcf_index}" ] && [ "~{sv_vcf_index}" != "" ]; then
        SV_IDX="~{sv_vcf_index}"
      fi
      if [ "~{sv_is_whole_genome}" = "true" ]; then
        prepare_vcf "~{sv_vcf}" sv.chr.vcf.gz true "${SV_IDX}"
      else
        prepare_vcf "~{sv_vcf}" sv.chr.vcf.gz false "${SV_IDX}"
      fi
      HAS_SV=1
    fi
    if [ -n "~{flare_vcf}" ] && [ "~{flare_vcf}" != "" ]; then
      prepare_vcf "~{flare_vcf}" flare.chr.vcf.gz false
      HAS_FLARE=1
    fi

    if [ "$((HAS_SNV + HAS_SV + HAS_FLARE))" -lt 1 ]; then
      echo "Need at least one of snv_indel_vcf, sv_vcf, flare_vcf" >&2
      exit 1
    fi

    HDR=""
    for f in snv_indel.chr.vcf.gz flare.chr.vcf.gz sv.chr.vcf.gz; do
      if [ -f "$f" ]; then HDR="$f"; break; fi
    done
    python3 - "$CHROM" "$HDR" <<'PY'
import re, sys, subprocess
chrom, hdr = sys.argv[1], sys.argv[2]
length = None
h = subprocess.check_output(["bcftools", "view", "-h", hdr], text=True)
pat = re.compile(r"^##contig=<ID=%s(?:,|>)" % re.escape(chrom))
for line in h.splitlines():
    if not pat.match(line):
        continue
    m = re.search(r"length=(\d+)", line)
    if m:
        length = int(m.group(1))
        break
FALLBACK = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
if length is None:
    length = FALLBACK.get(chrom) or FALLBACK.get(chrom.replace("chr", ""))
if length is None:
    raise SystemExit(f"could not determine length for {chrom}")
open("contig_length.txt", "w").write(str(length))
print(f"contig {chrom} length={length}", file=sys.stderr)
PY

    if [ "${HAS_SNV}" -eq 0 ]; then
      : > snv_indel.chr.vcf.gz
      : > snv_indel.chr.vcf.gz.tbi
    fi
    if [ "${HAS_SV}" -eq 0 ]; then
      : > sv.chr.vcf.gz
      : > sv.chr.vcf.gz.tbi
    fi
    if [ "${HAS_FLARE}" -eq 0 ]; then
      : > flare.chr.vcf.gz
      : > flare.chr.vcf.gz.tbi
    fi

    echo "${HAS_SNV}" | awk '{print ($1==1?"true":"false")}' > has_snv.txt
    echo "${HAS_SV}" | awk '{print ($1==1?"true":"false")}' > has_sv.txt
    echo "${HAS_FLARE}" | awk '{print ($1==1?"true":"false")}' > has_flare.txt
  >>>

  output {
    File snv_indel_chr_vcf = "snv_indel.chr.vcf.gz"
    File snv_indel_chr_tbi = "snv_indel.chr.vcf.gz.tbi"
    File sv_chr_vcf = "sv.chr.vcf.gz"
    File sv_chr_tbi = "sv.chr.vcf.gz.tbi"
    File flare_chr_vcf = "flare.chr.vcf.gz"
    File flare_chr_tbi = "flare.chr.vcf.gz.tbi"
    Int contig_length = read_int("contig_length.txt")
    Boolean has_snv = read_boolean("has_snv.txt")
    Boolean has_sv = read_boolean("has_sv.txt")
    Boolean has_flare = read_boolean("has_flare.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}

task MakeRegions {
  input {
    String chromosome
    Int contig_length
    Int n_shards
    # Prefer SNV/indel positions for density; fall back to FLARE then SV.
    File snv_indel_chr_vcf
    File snv_indel_chr_tbi
    File sv_chr_vcf
    File sv_chr_tbi
    File flare_chr_vcf
    File flare_chr_tbi
    Boolean has_snv
    Boolean has_sv
    Boolean has_flare
    String docker
    Int preemptible = 0
  }

  # Cromwell localizes all three chrom callsets; size disk from all of them.
  Int disk_gb = ceil(size(snv_indel_chr_vcf, "GB") + size(sv_chr_vcf, "GB") + size(flare_chr_vcf, "GB")) + 20

  command <<<
    set -euo pipefail
    CHROM="~{chromosome}"
    LENGTH=~{contig_length}
    N_SHARDS=~{n_shards}

    if [ "~{has_snv}" = "true" ]; then
      VCF="~{snv_indel_chr_vcf}"
      SRC=snv_indel
    elif [ "~{has_flare}" = "true" ]; then
      VCF="~{flare_chr_vcf}"
      SRC=flare
    elif [ "~{has_sv}" = "true" ]; then
      VCF="~{sv_chr_vcf}"
      SRC=sv
    else
      echo "Need at least one callset to build regions" >&2
      exit 1
    fi
    echo "MakeRegions density source=${SRC} file=$(basename "${VCF}")" >&2

    python3 - "$CHROM" "$LENGTH" "$N_SHARDS" "$VCF" <<'PY'
import subprocess
import sys

chrom, length_s, n_s, vcf = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
length = int(length_s)
n = max(1, int(n_s))

def n_records(path: str) -> int:
    try:
        out = subprocess.check_output(["bcftools", "index", "-n", path], text=True).strip()
        return int(out.split()[0])
    except Exception:
        # Fallback: stream-count (slow on huge files).
        p = subprocess.Popen(
            ["bcftools", "view", "-H", path],
            stdout=subprocess.PIPE,
            text=True,
        )
        assert p.stdout is not None
        c = sum(1 for _ in p.stdout)
        if p.wait() != 0:
            raise SystemExit(f"bcftools view failed on {path}")
        return c

total = n_records(vcf)
if total <= 0:
    # No sites — single full-contig shard.
    open("regions.txt", "w").write(f"{chrom}:1-{length}\n")
    print(f"wrote 1 region (0 sites in density VCF)", flush=True)
    raise SystemExit(0)

n = min(n, total)  # cannot have more shards than sites
# Balanced shard sizes (first total % n shards get one extra).
sizes = [total // n + (1 if i < (total % n) else 0) for i in range(n)]

proc = subprocess.Popen(
    ["bcftools", "query", "-f", "%POS\n", vcf],
    stdout=subprocess.PIPE,
    text=True,
)
assert proc.stdout is not None

ends = []  # inclusive end POS for shards 0..n-2
si = 0
remaining = sizes[0]
last_pos = 1
for line in proc.stdout:
    pos = int(line.strip())
    last_pos = pos
    remaining -= 1
    if remaining == 0 and si < n - 1:
        ends.append(pos)
        si += 1
        remaining = sizes[si]

rc = proc.wait()
if rc != 0:
    raise SystemExit(f"bcftools query failed (exit {rc})")

regions = []
start = 1
for end in ends:
    if end < start:
        # Pathological: should not happen with POS-sorted inputs.
        end = start
    regions.append(f"{chrom}:{start}-{end}")
    start = end + 1
regions.append(f"{chrom}:{start}-{length}")

# Drop inverted / empty windows (should be rare).
clean = []
for r in regions:
    _, span = r.split(":")
    a, b = map(int, span.split("-"))
    if a <= b:
        clean.append(r)
if not clean:
    clean = [f"{chrom}:1-{length}"]

open("regions.txt", "w").write("\n".join(clean) + "\n")
print(
    f"wrote {len(clean)} equal-variant regions for {chrom} "
    f"(sites={total}, sizes={sizes}, contig_len={length})",
    flush=True,
)
for r in clean:
    print(f"  {r}", flush=True)
PY
  >>>

  output {
    Array[String] regions = read_lines("regions.txt")
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}

task PropagateShard {
  input {
    String region
    File snv_indel_chr_vcf
    File snv_indel_chr_tbi
    File sv_chr_vcf
    File sv_chr_tbi
    File flare_chr_vcf
    File flare_chr_tbi
    Boolean has_snv
    Boolean has_sv
    Boolean has_flare
    String docker
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 3.0
    Int preemptible = 0
  }

  # Cromwell localizes full chrom callsets to every shard; we also write region
  # extracts + shard BCF. 3× covers inputs + extracts + output/headroom.
  Float input_gb = size(snv_indel_chr_vcf, "GB") + size(sv_chr_vcf, "GB") + size(flare_chr_vcf, "GB")
  Int disk_gb = ceil(input_gb * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    REGION="~{region}"
    echo "[$(date -Is)] PropagateShard ${REGION}" >&2
    df -h . || true

    subset() {
      local src="$1"
      local dst="$2"
      # Compressed VCF.gz — uncompressed -Ov with ~12k samples OOMs the disk.
      bcftools view -r "${REGION}" -Oz --threads ~{cpu} -o "${dst}" "${src}"
      bcftools index -t "${dst}"
      echo -n "$(basename "${dst}") records: "
      bcftools view -H "${dst}" | wc -l
      df -h . || true
    }

    ARGS=(
      --output shard.bcf
      --stats-json shard.stats.json
      --stats-tsv shard.stats.tsv
      --threads ~{cpu}
      --region "${REGION}"
      --progress-every 100000
    )

    if [ "~{has_snv}" = "true" ]; then
      subset "~{snv_indel_chr_vcf}" snv.region.vcf.gz
      ARGS+=(--snv-indel-vcf snv.region.vcf.gz)
    fi
    if [ "~{has_sv}" = "true" ]; then
      subset "~{sv_chr_vcf}" sv.region.vcf.gz
      ARGS+=(--sv-vcf sv.region.vcf.gz)
    fi
    if [ "~{has_flare}" = "true" ]; then
      subset "~{flare_chr_vcf}" flare.region.vcf.gz
      ARGS+=(--flare-vcf flare.region.vcf.gz)
    fi

    propagate-annotations "${ARGS[@]}"
    # Drop bulky region inputs before indexing output (peak disk is during merge).
    rm -f snv.region.vcf.gz snv.region.vcf.gz.tbi \
          sv.region.vcf.gz sv.region.vcf.gz.tbi \
          flare.region.vcf.gz flare.region.vcf.gz.tbi
    bcftools index -c shard.bcf
    df -h . || true
  >>>

  output {
    File shard_bcf = "shard.bcf"
    File shard_csi = "shard.bcf.csi"
    File stats_json = "shard.stats.json"
    File stats_tsv = "shard.stats.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}

task ConcatShards {
  input {
    Array[File] shard_bcfs
    Array[File] shard_stats_json
    String prefix
    Float? min_match_rate_flare = 0.95
    String docker
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 2.5
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(shard_bcfs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    echo "[$(date -Is)] ConcatShards" >&2

    cp ~{write_lines(shard_bcfs)} shard_list.txt
    bcftools concat -f shard_list.txt -Ob -o ~{prefix}.annotated.bcf --threads ~{cpu}
    bcftools index -c ~{prefix}.annotated.bcf

    export MIN_FLARE="~{min_match_rate_flare}"
    python3 - <<'PY'
import json, os, sys
from pathlib import Path
from collections import Counter

files = Path("~{write_lines(shard_stats_json)}").read_text().splitlines()
stats = [json.loads(Path(f).read_text()) for f in files if f.strip()]
if not stats:
    raise SystemExit("no shard stats")

def sum_key(key):
    return sum(s.get(key, 0) for s in stats)

merged = {
    "mode": "union_stream_sharded",
    "output_vcf": "~{prefix}.annotated.bcf",
    "chromosome": stats[0].get("chromosome", ""),
    "n_samples": stats[0].get("n_samples", 0),
    "union_sites": sum_key("union_sites"),
    "output_sites": sum_key("output_sites"),
    "sites_by_n_sources": {},
    "snv_flare_overlap_rate": None,
    "sources": [],
    "n_shards": len(stats),
}

hist = Counter()
for s in stats:
    for k, v in (s.get("sites_by_n_sources") or {}).items():
        hist[str(k)] += v
merged["sites_by_n_sources"] = {int(k): v for k, v in sorted(hist.items(), key=lambda kv: int(kv[0]))}

labels = [src["label"] for src in stats[0].get("sources", [])]
for label in labels:
    pieces = []
    for s in stats:
        for src in s.get("sources", []):
            if src["label"] == label:
                pieces.append(src)
                break
    source_sites = sum(p.get("source_sites", 0) for p in pieces)
    with_src = sum(p.get("union_sites_with_source", 0) for p in pieces)
    gt_filled = sum(p.get("sample_gt_filled", 0) for p in pieces)
    gt_total = sum(p.get("sample_gt_total", 0) for p in pieces)
    union = merged["union_sites"]
    merged["sources"].append({
        "label": label,
        "path": pieces[0].get("path", "") if pieces else "",
        "source_sites": source_sites,
        "union_sites_with_source": with_src,
        "union_sites_without_source": max(0, union - with_src),
        "union_annotation_rate": (with_src / union) if union else 0.0,
        "sample_gt_fill_rate": (gt_filled / gt_total) if gt_total else 0.0,
        "sample_gt_filled": gt_filled,
        "sample_gt_total": gt_total,
    })

snv_sites = sum(s.get("snv_sites", 0) for s in stats)
flare_on_snv = sum(s.get("flare_on_snv", 0) for s in stats)
if snv_sites:
    merged["snv_flare_overlap_rate"] = flare_on_snv / snv_sites
merged["snv_sites"] = snv_sites
merged["flare_on_snv"] = flare_on_snv

prefix = "~{prefix}"
Path(f"{prefix}.propagate_stats.json").write_text(json.dumps(merged, indent=2) + "\n")

lines = ["metric\tvalue",
         f"mode\t{merged['mode']}",
         f"chromosome\t{merged['chromosome']}",
         f"n_samples\t{merged['n_samples']}",
         f"n_shards\t{merged['n_shards']}",
         f"union_sites\t{merged['union_sites']}",
         f"output_sites\t{merged['output_sites']}"]
if merged["snv_flare_overlap_rate"] is not None:
    lines.append(f"snv_flare_overlap_rate\t{merged['snv_flare_overlap_rate']:.6f}")
for n, c in merged["sites_by_n_sources"].items():
    lines.append(f"union_sites_with_{n}_sources\t{c}")
for src in merged["sources"]:
    p = src["label"]
    lines += [
        f"{p}_source_sites\t{src['source_sites']}",
        f"{p}_union_sites_with_source\t{src['union_sites_with_source']}",
        f"{p}_union_sites_without_source\t{src['union_sites_without_source']}",
        f"{p}_union_annotation_rate\t{src['union_annotation_rate']:.6f}",
        f"{p}_sample_gt_fill_rate\t{src['sample_gt_fill_rate']:.6f}",
    ]
Path(f"{prefix}.propagate_stats.tsv").write_text("\n".join(lines) + "\n")

threshold = os.environ.get("MIN_FLARE", "").strip()
if threshold and threshold.lower() not in ("", "null") and merged["snv_flare_overlap_rate"] is not None:
    thr = float(threshold)
    rate = merged["snv_flare_overlap_rate"]
    if rate < thr:
        raise SystemExit(f"flare_on_snv match rate too low: {rate:.4f} < {thr:.4f}")

print("=== merged propagate stats ===")
print(Path(f"{prefix}.propagate_stats.tsv").read_text())
PY
    echo -n "annotated records: "
    bcftools view -H ~{prefix}.annotated.bcf | wc -l
    echo "[$(date -Is)] ConcatShards done" >&2
  >>>

  output {
    File annotated_bcf = "~{prefix}.annotated.bcf"
    File annotated_bcf_csi = "~{prefix}.annotated.bcf.csi"
    File propagate_stats_json = "~{prefix}.propagate_stats.json"
    File propagate_stats_tsv = "~{prefix}.propagate_stats.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}

workflow PropagateAnnotations {
  input {
    String chromosome
    File? snv_indel_vcf
    File? sv_vcf
    File? sv_vcf_index
    File? flare_vcf
    String prefix = "annotated"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-propagate-annotations:0.3.1"
    Float? min_match_rate_flare = 0.95
    Boolean sv_is_whole_genome = true
    Int n_shards = 8
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 2.0
    Int shard_disk_gb_floor = 50
    Float shard_disk_gb_multiplier = 3.0
    Int preemptible = 0
  }

  call PrepareChrom {
    input:
      chromosome = chromosome,
      snv_indel_vcf = snv_indel_vcf,
      sv_vcf = sv_vcf,
      sv_vcf_index = sv_vcf_index,
      flare_vcf = flare_vcf,
      sv_is_whole_genome = sv_is_whole_genome,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier,
      preemptible = preemptible
  }

  call MakeRegions {
    input:
      chromosome = chromosome,
      contig_length = PrepareChrom.contig_length,
      n_shards = n_shards,
      snv_indel_chr_vcf = PrepareChrom.snv_indel_chr_vcf,
      snv_indel_chr_tbi = PrepareChrom.snv_indel_chr_tbi,
      sv_chr_vcf = PrepareChrom.sv_chr_vcf,
      sv_chr_tbi = PrepareChrom.sv_chr_tbi,
      flare_chr_vcf = PrepareChrom.flare_chr_vcf,
      flare_chr_tbi = PrepareChrom.flare_chr_tbi,
      has_snv = PrepareChrom.has_snv,
      has_sv = PrepareChrom.has_sv,
      has_flare = PrepareChrom.has_flare,
      docker = docker,
      preemptible = preemptible
  }

  scatter (region in MakeRegions.regions) {
    call PropagateShard {
      input:
        region = region,
        snv_indel_chr_vcf = PrepareChrom.snv_indel_chr_vcf,
        snv_indel_chr_tbi = PrepareChrom.snv_indel_chr_tbi,
        sv_chr_vcf = PrepareChrom.sv_chr_vcf,
        sv_chr_tbi = PrepareChrom.sv_chr_tbi,
        flare_chr_vcf = PrepareChrom.flare_chr_vcf,
        flare_chr_tbi = PrepareChrom.flare_chr_tbi,
        has_snv = PrepareChrom.has_snv,
        has_sv = PrepareChrom.has_sv,
        has_flare = PrepareChrom.has_flare,
        docker = docker,
        cpu = cpu,
        memory_gb = memory_gb,
        disk_gb_floor = shard_disk_gb_floor,
        disk_gb_multiplier = shard_disk_gb_multiplier,
        preemptible = preemptible
    }
  }

  call ConcatShards {
    input:
      shard_bcfs = PropagateShard.shard_bcf,
      shard_stats_json = PropagateShard.stats_json,
      prefix = prefix,
      min_match_rate_flare = min_match_rate_flare,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier,
      preemptible = preemptible
  }

  output {
    File annotated_bcf = ConcatShards.annotated_bcf
    File annotated_bcf_csi = ConcatShards.annotated_bcf_csi
    File propagate_stats_json = ConcatShards.propagate_stats_json
    File propagate_stats_tsv = ConcatShards.propagate_stats_tsv
    Array[File] shard_bcfs = PropagateShard.shard_bcf
    Int n_regions = length(MakeRegions.regions)
  }

  meta {
    description: "Sharded per-chromosome union of SNV/indel, SV, and FLARE into one annotated BCF."
  }
}
