# Phase 3 exact panel sampling-size sensitivity

## Purpose

The replicated `1 Mb` exact panel pilot showed:

- stable `Campodea` preference for symmetric metrics;
- unstable directional `containment_a_in_b` and overlap-style behavior for `protura_purged`.

That left an open question:

- was the directional instability mostly caused by insufficient sampled sequence?

This follow-up pilot repeats the same deterministic window-sampling design at `5 Mb` per taxon to test whether larger sampled sequence stabilizes the rankings.

This remains a methodological size-sensitivity pilot only.

## Design

Taxa:

- `protura_purged`
- `GCA_009757345.1` / *Campodea augens*
- `GCA_000934665.2` / *Catajapyx aquilonaris*
- `GCA_052721305.1` / *Folsomia candida*
- `GCA_020796165.1` / *Pantala flavescens*
- `GCA_021134715.1` / *Daphnia pulex*

Source FASTA type:

- masked Phase 2 FASTA
- `results/phase2_expanded_panel/repeats/masked/<sample>.fa`

Sampling:

- same deterministic window-based method as the `1 Mb` pilot
- non-overlapping `50 kb` windows
- fixed seed per replicate and sample
- windows concatenated until `5,000,000 bp`
- only the final selected window truncated if needed

Output root:

- `results/phase3_metric_panel_pilot/replicates_5mb/`

Replicates:

- `5`

k-mer settings:

- `k = 15`
- `canonical = true`

## Completion and practical notes

All `30` taxon-replicate `5 Mb` mini-FASTA subsets completed successfully.

Manifest:

- `results/phase3_metric_panel_pilot/replicates_5mb/mini_fasta_sampling_manifest.tsv`

Per-replicate exact panel outputs:

- `results/phase3_metric_panel_pilot/replicates_5mb/exact_panel_k15_5mb_rep01.tsv`
- `results/phase3_metric_panel_pilot/replicates_5mb/exact_panel_k15_5mb_rep02.tsv`
- `results/phase3_metric_panel_pilot/replicates_5mb/exact_panel_k15_5mb_rep03.tsv`
- `results/phase3_metric_panel_pilot/replicates_5mb/exact_panel_k15_5mb_rep04.tsv`
- `results/phase3_metric_panel_pilot/replicates_5mb/exact_panel_k15_5mb_rep05.tsv`

Summary outputs:

- `results/phase3_metric_panel_pilot/replicates_5mb/replicate_nearest_by_metric_k15_5mb.tsv`
- `results/phase3_metric_panel_pilot/replicates_5mb/aggregate_nearest_by_metric_k15_5mb.tsv`

Runtime note:

- exact `5 Mb` panel runs were noticeably slower than the `1 Mb` pilot, as expected for full in-memory unique-kmer sets
- this remains practical for small pilots, but it reinforces that the exact-panel script is still a validation tool rather than a scalable production path

## protura_purged at 5 Mb

Aggregate nearest-neighbor pattern for `protura_purged`:

| metric | top match pattern across 5 replicates |
| --- | --- |
| `jaccard` | `Campodea` in `5/5` |
| `dice` | `Campodea` in `5/5` |
| `binary_cosine` | `Campodea` in `4/5`, `Pantala` in `1/5` |
| `mash_distance` | `Campodea` in `5/5` |
| `containment_b_in_a` | `Campodea` in `5/5` |
| `containment_a_in_b` | `Pantala` in `4/5`, `Folsomia` in `1/5` |
| `overlap_coefficient` | `Pantala` in `4/5`, `Folsomia` in `1/5` |

Representative `5 Mb` aggregate rows for `protura_purged`:

- `jaccard -> Campodea`: `5/5`, mean `0.014632539663`
- `dice -> Campodea`: `5/5`, mean `0.028842973430`
- `binary_cosine -> Campodea`: `4/5`, mean `0.028897321719`
- `mash_distance -> Campodea`: `5/5`, mean `0.236397259210`
- `containment_b_in_a -> Campodea`: `5/5`, mean `0.028322244254`
- `containment_a_in_b -> Pantala`: `4/5`, mean `0.031569846046`
- `overlap_coefficient -> Pantala`: `4/5`, mean `0.031569846046`

## Comparison with 1 Mb

### Symmetric metrics

At `1 Mb`:

- `jaccard`, `dice`, `binary_cosine`, and `mash_distance` all preferred `Campodea` in `5/5`

At `5 Mb`:

- `jaccard`, `dice`, and `mash_distance` still preferred `Campodea` in `5/5`
- `binary_cosine` remained almost fully stable, with `Campodea` in `4/5` and `Pantala` in `1/5`

Conclusion:

- the symmetric metrics remain strongly stable with more sampled sequence
- larger sampled sequence mainly increases absolute similarity values, not the qualitative ranking pattern

### containment_b_in_a

At `1 Mb`:

- `Campodea` in `5/5`

At `5 Mb`:

- `Campodea` in `5/5`

Conclusion:

- the directional “reference content recovered in protura_purged” direction remains fully stable and aligned with the symmetric metrics

### containment_a_in_b

At `1 Mb`:

- `Pantala` in `3/5`
- `Campodea` in `1/5`
- `Folsomia` in `1/5`

At `5 Mb`:

- `Pantala` in `4/5`
- `Folsomia` in `1/5`
- `Campodea` no longer appears as the top match

Conclusion:

- larger sampled sequence does not stabilize `containment_a_in_b` toward the symmetric `Campodea` pattern
- instead, it makes `Pantala` more consistently preferred

### overlap_coefficient

At `1 Mb`:

- `Campodea` in `2/5`
- `Pantala` in `2/5`
- `Folsomia` in `1/5`

At `5 Mb`:

- `Pantala` in `4/5`
- `Folsomia` in `1/5`
- `Campodea` disappears as the top overlap match

Conclusion:

- larger sampled sequence does not collapse overlap behavior toward the symmetric ranking
- it pushes overlap toward the same side as `containment_a_in_b`

## Interpretation

### Does larger sampled sequence stabilize containment/overlap?

Only partially, and not in the direction that would reduce the metric disagreement.

Observed effect:

- symmetric metrics were already stable at `1 Mb` and remain stable at `5 Mb`
- `containment_a_in_b` and overlap become more internally consistent at `5 Mb`
- but they stabilize around `Pantala` rather than converging toward `Campodea`

This means the earlier directional instability was not just a low-coverage artifact of the `1 Mb` subset size.

### Does metric choice still matter at 5 Mb?

Yes, clearly.

At `5 Mb`:

- symmetric metrics and `containment_b_in_a` still favor `Campodea`
- `containment_a_in_b` and overlap now favor `Pantala` even more consistently

So the disagreement between metric families persists after increasing sampled sequence fivefold.

## Conservative conclusion

The `5 Mb` pilot strengthens the methodological conclusion:

- metric choice, not just subset size, is a major driver of ranking behavior

What changes from `1 Mb` to `5 Mb`:

- absolute scores increase substantially
- directional/smaller-set-sensitive metrics become less noisy

What does not change:

- the core split between symmetric metrics and `containment_a_in_b` / overlap remains

## Does this support implementing containment in production now?

It supports continued interest in containment-like metrics, but still argues for caution.

Why caution is still warranted:

- directional containment is clearly measuring something different from symmetric overlap
- that difference persists even after more sequence is sampled
- the result is now more stable, but it still does not tell us which quantity best matches the intended biological or query/reference use case

Conservative implication:

- this is now strong evidence that containment-like metrics deserve explicit evaluation as a separate production option
- it is not yet sufficient evidence to replace or reinterpret the current symmetric sketch branch without a clearer statement of the target use case

## Limitations

- only six taxa
- only masked FASTA
- only `k = 15`
- only `1 Mb` and `5 Mb` subset sizes
- exact sets remain in-memory and non-scalable for large panels
- no trees were built
- no whole-genome exact comparisons were attempted
- no claim is made about biological truth
