# The Mendelian violation rate is not comparable across `nanc`

Why `sel_chr20_flare2_nanc5` and `sel_chr20_flare2_nanc6` cannot be ranked by
the raw selection-v2 metric, and what to do instead.

Written 2026-10-08. Nothing here uses AoU data; it is a property of the metric.

---

## 1. Two separate problems

### 1.1 A silent bug: codes above 4 were read as missing

`scripts/flare_score_mendelian_lai.py` hard-coded
`ANCESTRY = {0: eas, 1: amr, 2: eur, 3: afr, 4: sas}` and dropped any ancestry
call not in it. A FLARE2 model with `nanc` 6 emits code 5. A locus is scored
only when all six haplotype calls are present, so **every locus touching
cluster 5 fell out of `n_informative_locus_calls`** — no warning, no counter.

Measured on a synthetic trio with `nanc` 6 and uniform cluster usage: 400 loci
in, **196 informative**, exactly `400 − 204` where 204 is the number of loci
touching code 5. Over half the chromosome, discarded non-randomly — and
precisely at the loci where the sixth cluster is doing work.

This is the same failure mode the 2026-10-02 review already caught once
("the Mendelian gate never ran … the code skipped silently").

**Fixed.** The scorer now reads the VCF's `##ANCESTRY=<…>` header (the
convention `flare_score_allele_ancestry.py` already follows), accepts
`--num-ancs`, counts out-of-range calls, and **refuses to emit a score** when
it sees any, unless `--allow-out-of-range` is passed. The summary records
`ancestry_alphabet`, `ancestry_alphabet_source`, `n_ancestries` and
`grid_call_stats`. Regression tests:
`scripts/test_flare_score_mendelian_nanc.py`.

### 1.2 A structural bias: the metric rewards fewer labels

A locus is a hard violation when the child's two haplotype labels cannot be
drawn one from each parent. **Compatibility is preserved under coarsening**: if
two labels are merged, every compatible locus stays compatible and some
violations disappear. Nothing works the other way. So at fixed painting
accuracy the violation rate rises with the number of labels, independent of
whether the painting is better.

`simulate_k_bias.py` quantifies it using the production `score_trio_path`, on
101 trios (the Phase 2 pedigree), 3000 loci, 3 replicates:

| LAI error per haplotype-locus | K=3 | K=5 | K=6 | K=8 |
|---|---|---|---|---|
| 0% | 0 | 0 | 0 | 0 |
| 1% | 2.17% | 2.91% | 3.14% | 3.33% |
| 2% | 4.28% | 5.71% | 6.14% | 6.55% |
| 5% | 10.2% | 13.8% | 14.7% | 15.7% |
| 10% | 18.9% | 25.6% | 27.3% | 29.3% |

K=6 is **7.7% worse than K=5 at identical accuracy** (2% error). That alone is
comparable to the differences a recipe comparison is trying to resolve.

The sharper test holds the truth fixed at 5 ancestries and lets the K=6 model
add one noisy label — a cluster that captures no new structure, which is what
the chr20 FLARE2 fits may well be doing (`flare/README.md`: the sixth cluster
is SAS-like with EAS/EUR/SAS weights near 1/3 each, and both fits carry a
low-autocorrelation cluster at ~1.5–1.7% in every population):

| Split purity | K=5 | K=6 as-is | K=6 projected to panels |
|---|---|---|---|
| 85% | 5.4% | **18.1%** | 4.8% |
| 95% | 5.3% | **10.4%** | 4.8% |

A label carrying no information inflates the rate **2–3.4×**. Projecting it
back onto its dominant panel removes the inflation.

---

## 2. What this means for the nanc 5 vs 6 decision

1. **The raw ranking cannot decide it.** Both effects push the same way, and
   both favour `nanc5` — and also favour the K=5 original-FLARE pin
   (`chr20_flat_props_pin_t`) over `nanc6`. A raw-rate win for `nanc5` or the
   pin is uninformative about `nanc6`. A raw-rate win for **`nanc6`** would be
   meaningful, because it would have won against a headwind.
2. **Score every candidate at one resolution.** Pass `--project-labels` so each
   cluster is replaced by its `dominant_panel`
   (`flare2_build_model.py` already writes this to `<prefix>.labels.tsv`), and
   all candidates land in the five-panel alphabet. For `nanc` 6 the local
   `nanc6.model` relabels to `afr, amr, eas, eas, eur, eur`, i.e. four distinct
   panels; for the original-FLARE pin the ancestries are the panels already.
3. **`select_recipe` now refuses the mistake.** It compares `projected_to` /
   `ancestry_alphabet` across survivors and returns
   `metric_incomparable_mixed_ancestry_alphabets` with no winner rather than
   ranking across resolutions.

The projection is not free: merging two clusters can mask a real inconsistency
between them, and in the split simulation it makes K=6 look marginally *better*
than K=5 (4.8% vs 5.3–5.4%), because an error landing on the sibling label
becomes invisible. That residual is ~10% relative, against the 2–3.4× it
removes. Treat the projected rate as the ranking metric and the raw rate as a
within-recipe diagnostic only.

---

## 3. Caveats

- The simulation is a Markov admixture mosaic (switch rate `T`=10 per Morgan,
  one crossover per Morgan, uniform misclassification). It is calibrated to the
  real trio count, not to AoU painting behaviour; the numbers bound the
  direction and rough size of the bias, they are not predictions of what chr20
  will give.
- `nanc5` and `nanc6` are independent GMM fits, not an exact
  coarsening/refinement pair, so the monotonicity argument in §1.2 is a
  mechanism, not a theorem about these two specific models.
- The misclassification model is uniform-over-other-labels. Real LAI errors
  concentrate between similar ancestries, which would reduce the bias somewhat
  for clusters that are genuinely distinct and increase it for clusters that
  are not.
- **The size of the effect depends on label frequencies, not only on K.** The
  simulation above uses uniform label frequencies (+9-11% for K=4 to 5). Scoring
  one synthetic pin-like VCF with and without folding a rare label (SAS, 5% of
  ancestry) into EUR gave a five-label rate 23-26% higher
  (`fold_rare_label.py`). Part 8 uses the larger figure as the tie allowance
  (`LABEL_COUNT_BIAS_PER_LABEL = 0.26` in `scripts/flare_lai_exp.py`) when it
  ranks the five-label pin against the four-label FLARE2 recipes, and marks a
  lead taken by the fewer-label recipe provisional. Scoring both on one alphabet
  (`--project-labels`, `PIN_SAS_TO` in cell 26) is the only exact remedy.

---

## 4. Rerun

```bash
python3 flare/eval/mendelian_k_bias/simulate_k_bias.py --out-dir flare/eval/mendelian_k_bias
python3 flare/eval/mendelian_k_bias/plot_k_bias.py \
  --in-dir flare/eval/mendelian_k_bias --out flare/eval/mendelian_k_bias/mendelian_k_bias.png
python3 flare/eval/mendelian_k_bias/fold_rare_label.py
python3 -m pytest scripts/test_flare_score_mendelian_nanc.py scripts/test_flare_selection_v2.py \
  scripts/test_flare_alphabet_pipeline.py -q
```

Outputs: `k_bias_raw.tsv`, `k_bias_summary.tsv`, `k_bias_summary.json`,
`mendelian_k_bias.png`.
