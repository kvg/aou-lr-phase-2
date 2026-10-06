# Figure 2 variant length spectrum

1 bp signed-length counts, split by sequence context (US / RM / SD / SR), built
from the same sources and rules as Table 2 (`tab:callset`). Display binning and
facet cuts are applied later when rendering Figure 2.

## Sources

| Partition | Source | Rule (matches Table 2) |
|---|---|---|
| `small` | DeepVariant + GLnexus shards (`GL_INTERVAL_set`, `this.VCF`) | FILTER `PASS` or `.`; SNVs (length 0) and indels with \|L\| < 20 |
| `sv` | `v3_main.bcf` | Resolved DEL / INS with \|SVLEN\| ≥ 20 (`size_bin_ge20`) |
| `ultralong` | `v3_ultralong.bcf` | Every record |
| BND | `v3_bnd.bcf` | Count only (no length) |

The `small` → `sv` source switch at 20 bp is the one Table 2 already uses.
GLnexus indels ≥ 20 bp and main SVs < 20 bp are counted in the summary
(`glnexus_indel_ge_small_max_dropped`, `main_lt_sv_min_dropped`) but not in
the rows. MNPs, spanning `*` and symbolic GLnexus alleles are summary only.
SVTYPE / SVLEN / END for the SV companions are derived exactly as in
`extract_sites_from_vcf.py`; companions are counted genome-wide, as in Table 2.

`region_class` uses `annotate_repeat_context.context_label` with the merged
`hg38.RM/SR/SD` tracks: either breakpoint assigns the class (SR > SD > RM),
DEL/DUP/CNV > 5 kb use body coverage > 0.5, insertions collapse to the anchor.
Small variants use the same breakpoint rule on their REF span. Use the same
BED URIs as the AnnotateSvCallset run behind Table 2.

## Outputs

| File | Contents |
|---|---|
| `*.length_spectrum.bins.tsv` | `chrom partition signed_len region_class n_sites` (sparse 1 bp) |
| `*.length_spectrum.summary.json` | `counts` (record / allele tallies); merged summary adds `table2_check` |

## Terra launch

1. Sync scripts + WDL: `python3 scripts/terra_sync_repo.py`.
2. [`wdl/VariantLengthSpectrum.wdl`](wdl/VariantLengthSpectrum.wdl), root entity
   `GL_INTERVAL_set`, one row per chromosome (chr1–22, X, Y; skip chrM).
   Inputs: [`configs/variant_length_spectrum.inputs.json.example`](configs/variant_length_spectrum.inputs.json.example).
3. [`wdl/MergeVariantLengthSpectrum.wdl`](wdl/MergeVariantLengthSpectrum.wdl)
   on the collected shard outputs. This also scans `v3_main` / ultralong / BND
   once. Inputs: [`configs/merge_variant_length_spectrum.inputs.json.example`](configs/merge_variant_length_spectrum.inputs.json.example)
   (or pass the shard output arrays directly).

Image: `aou-sv-annotation:0.1.6` (bcftools). Scripts are localized as task
inputs, so no Docker rebuild is needed.

Check `table2_check` in the merged summary against Table 2 before plotting,
then copy the merged TSV / JSON to
`aou-lr-phase-2-manuscript/mockup_figures/data/fig2_length_spectrum.*`.
