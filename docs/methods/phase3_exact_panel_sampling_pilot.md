# Phase 3 exact panel sampling pilot

## Purpose

The first exact mini-FASTA pilot used only the first `1 Mb` per taxon. That was enough to show a metric-dependent ranking difference for `protura_purged`, but not enough to tell whether the difference reflected:

- the metric itself;
- or
- arbitrary sequence-subset choice.

This follow-up pilot therefore adds replicated mini-FASTA sampling on the same six taxa to test ranking stability under repeated sequence subset selection.

This remains a methodological sampling validation only.

## Taxa and source FASTA type

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

## Sampling design

Replicates:

- `5`

Cap per replicate:

- `1,000,000 bp` per taxon

k-mer settings:

- `k = 15`
- `canonical = true`

Sampling method:

- all masked FASTA contigs were partitioned into non-overlapping windows of `50 kb`
- final remainder windows shorter than `50 kb` were retained
- for each taxon and replicate, the window list was deterministically shuffled with a fixed seed:
  - seed template: `2026-07-01-phase3b:<replicate>:<sample>`
- shuffled windows were concatenated until `1 Mb` was reached
- if the last selected window exceeded the remaining cap, that final window was truncated cleanly

This means the pilot samples windows rather than full contigs.

Why windows were used:

- they avoid always taking the first contigs;
- they allow reproducible replicated subsets even for taxa dominated by a few large contigs;
- they still preserve real sequence blocks and contig boundaries at the window level.

Manifest:

- `results/phase3_metric_panel_pilot/replicates/mini_fasta_sampling_manifest.tsv`

Observed extraction summary:

- all `30` taxon-replicate mini-FASTA files completed successfully
- every taxon-replicate reached exactly `1,000,000 bp`
- `Pantala`, `Daphnia`, and `Folsomia` typically used `20` full windows with no truncation
- `protura_purged`, `Campodea`, and `Catajapyx` often required more windows and a truncated final window

## Exact panel runs

Per replicate, the standalone exact panel was run on the six mini-FASTA files:

- `results/phase3_metric_panel_pilot/replicates/exact_panel_k15_1mb_rep01.tsv`
- `results/phase3_metric_panel_pilot/replicates/exact_panel_k15_1mb_rep02.tsv`
- `results/phase3_metric_panel_pilot/replicates/exact_panel_k15_1mb_rep03.tsv`
- `results/phase3_metric_panel_pilot/replicates/exact_panel_k15_1mb_rep04.tsv`
- `results/phase3_metric_panel_pilot/replicates/exact_panel_k15_1mb_rep05.tsv`

Summary outputs:

- `results/phase3_metric_panel_pilot/replicates/replicate_nearest_by_metric_k15_1mb.tsv`
- `results/phase3_metric_panel_pilot/replicates/aggregate_nearest_by_metric_k15_1mb.tsv`

## protura_purged summary across replicates

Nearest-neighbor stability for `protura_purged`:

| metric | top match pattern across 5 replicates |
| --- | --- |
| `jaccard` | `Campodea` in `5/5` |
| `dice` | `Campodea` in `5/5` |
| `binary_cosine` | `Campodea` in `5/5` |
| `mash_distance` | `Campodea` in `5/5` |
| `containment_b_in_a` | `Campodea` in `5/5` |
| `containment_a_in_b` | `Pantala` in `3/5`, `Campodea` in `1/5`, `Folsomia` in `1/5` |
| `overlap_coefficient` | `Campodea` in `2/5`, `Pantala` in `2/5`, `Folsomia` in `1/5` |

Representative aggregate rows for `protura_purged`:

- `jaccard -> Campodea`: `5/5`, mean `0.003502567251`
- `dice -> Campodea`: `5/5`, mean `0.006980648177`
- `binary_cosine -> Campodea`: `5/5`, mean `0.006980864208`
- `mash_distance -> Campodea`: `5/5`, mean `0.331022794098`
- `containment_b_in_a -> Campodea`: `5/5`, mean `0.006997507404`
- `containment_a_in_b -> Pantala`: `3/5`, mean `0.007355102126`
- `overlap_coefficient -> mixed`: `Campodea 2/5`, `Pantala 2/5`, `Folsomia 1/5`

## Interpretation

### Was Campodea versus Pantala stable?

For `protura_purged`:

- yes for the symmetric metrics tested here:
  - Jaccard
  - Dice
  - binary cosine
  - Mash distance
- yes for directional `containment_b_in_a`
- no for directional `containment_a_in_b`
- no for overlap coefficient

This means the first pilot’s `Campodea` signal for symmetric metrics was not just an artifact of taking the first `1 Mb`.

### Does metric choice still matter after resampling?

Yes.

That conclusion is stronger after replication than before:

- symmetric metrics remained consistently `Campodea`-favoring for `protura_purged`
- directional containment from the `protura_purged` side remained unstable and often favored `Pantala`
- overlap coefficient also showed sampling sensitivity and did not settle on a single stable top match

### Is instability mostly metric or sampling?

For `protura_purged`, the dominant effect appears to be metric choice, with sampling as a secondary effect.

Reasoning:

- if sampling were the main driver, symmetric metrics would also be expected to fluctuate more strongly across replicates
- instead, the symmetric metrics were fully stable across `5/5` replicates
- the instability is concentrated in:
  - `containment_a_in_b`
  - overlap coefficient

That pattern is more consistent with different metrics emphasizing different aspects of overlap, with low absolute scores making some directional/smaller-set style metrics more sensitive to subset choice.

## Conservative conclusion

This replicated pilot supports two points:

1. the metric-dependent disagreement observed in the first mini-FASTA pilot is real enough to persist after resampling
2. some metrics, especially directional containment and overlap-style quantities, are themselves sensitive to which real sequence subset is sampled

This does not justify biological interpretation.

It does justify continued methodological evaluation before changing the production sketch branch.

## Limitations

- only six taxa were used
- only masked FASTA was tested
- only `1 Mb` per taxon was sampled
- windows, not whole contigs or whole assemblies, were sampled
- exact in-memory sets were used, so this is still not a scalable whole-panel method
- no trees were built
- no whole-genome exact comparison was attempted
- absolute similarity values remained very low, so ranking differences should be treated cautiously

## Practical implication

At this stage, a production change should not be justified by one pilot ranking alone.

The evidence is now strong enough to say:

- symmetric exact overlap metrics and directional containment are not interchangeable
- their behavior on real data can remain different even after replicated resampling
