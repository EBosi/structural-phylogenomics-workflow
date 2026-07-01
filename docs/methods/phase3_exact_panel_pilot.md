# Phase 3B exact metric panel pilot

## Purpose

This pilot uses the standalone exact k-mer distance panel on tiny real Phase 2 FASTA extracts to test whether different set-based metrics produce different nearest-neighbor behavior on real sequence content.

This is a methodological pilot only.

It is not:

- a whole-genome analysis;
- a replacement for the Phase 2 sketch workflow;
- a biological conclusion;
- a tree-based comparison.

## Input subset design

Taxa included:

- `protura_purged`
- `GCA_009757345.1` / *Campodea augens*
- `GCA_000934665.2` / *Catajapyx aquilonaris*
- `GCA_052721305.1` / *Folsomia candida*
- `GCA_020796165.1` / *Pantala flavescens*
- `GCA_021134715.1` / *Daphnia pulex*

Source FASTA type used:

- masked Phase 2 FASTA
- specifically `results/phase2_expanded_panel/repeats/masked/<sample>.fa`

Reason for using masked:

- masked high-k inputs were already produced uniformly for these taxa;
- this avoids rerunning preprocessing and keeps the pilot aligned with the existing high-k branch input family.

## Source FASTA paths used

- `results/phase2_expanded_panel/repeats/masked/protura_purged.fa`
- `results/phase2_expanded_panel/repeats/masked/GCA_009757345.1.fa`
- `results/phase2_expanded_panel/repeats/masked/GCA_000934665.2.fa`
- `results/phase2_expanded_panel/repeats/masked/GCA_052721305.1.fa`
- `results/phase2_expanded_panel/repeats/masked/GCA_020796165.1.fa`
- `results/phase2_expanded_panel/repeats/masked/GCA_021134715.1.fa`

## Mini-FASTA creation

Pilot output directory:

- `results/phase3_metric_panel_pilot/`

Mini-FASTA directory:

- `results/phase3_metric_panel_pilot/mini_fastas/`

Cap size:

- `1,000,000 bp` total per taxon

Mini-FASTA extraction rule:

- preserve contig boundaries while possible;
- append full contigs until the remaining budget is smaller than the next contig;
- truncate the last included contig cleanly to fill the remaining cap;
- stop at exactly `1,000,000 bp`.

Manifest:

- `results/phase3_metric_panel_pilot/mini_fasta_manifest.tsv`

Observed extraction summary:

- all six taxa reached the `1 Mb` cap;
- all six mini-FASTA files required truncation of the last included contig;
- contig counts written ranged from `1` to `5`.

## Command used

```bash
python workflow/scripts/compute_kmer_set_distance_panel.py \
  --k 15 \
  --canonical true \
  --output results/phase3_metric_panel_pilot/exact_panel_k15_1mb.tsv \
  results/phase3_metric_panel_pilot/mini_fastas/protura_purged.fa \
  results/phase3_metric_panel_pilot/mini_fastas/GCA_009757345.1.fa \
  results/phase3_metric_panel_pilot/mini_fastas/GCA_000934665.2.fa \
  results/phase3_metric_panel_pilot/mini_fastas/GCA_052721305.1.fa \
  results/phase3_metric_panel_pilot/mini_fastas/GCA_020796165.1.fa \
  results/phase3_metric_panel_pilot/mini_fastas/GCA_021134715.1.fa
```

Main outputs:

- `results/phase3_metric_panel_pilot/exact_panel_k15_1mb.tsv`
- `results/phase3_metric_panel_pilot/nearest_by_metric_k15_1mb.tsv`

## Nearest-neighbor summary

For `protura_purged`, the best non-self match by metric was:

| metric | best match | score |
| --- | --- | ---: |
| `jaccard` | `GCA_009757345.1` / *Campodea augens* | `0.0036370390` |
| `dice` | `GCA_009757345.1` / *Campodea augens* | `0.0072477178` |
| `overlap_coefficient` | `GCA_009757345.1` / *Campodea augens* | `0.0076798648` |
| `containment_a_in_b` | `GCA_020796165.1` / *Pantala flavescens* | `0.0075174845` |
| `containment_b_in_a` | `GCA_009757345.1` / *Campodea augens* | `0.0076798648` |
| `binary_cosine` | `GCA_009757345.1` / *Campodea augens* | `0.0072592195` |
| `mash_distance` | `GCA_009757345.1` / *Campodea augens* | `0.3284712430` |

Ranking for `protura_purged`:

- Jaccard-like symmetric metrics ranked:
  1. `Campodea`
  2. `Pantala`
  3. `Folsomia`
  4. `Daphnia`
  5. `Catajapyx`
- directional `containment_a_in_b` ranked:
  1. `Pantala`
  2. `Campodea`
  3. `Folsomia`
  4. `Daphnia`
  5. `Catajapyx`

## Did Jaccard and containment disagree?

Yes, but only modestly in this mini-FASTA pilot.

Observed for `protura_purged`:

- `jaccard`, `dice`, `overlap_coefficient`, `binary_cosine`, and `mash_distance` all preferred `Campodea`
- `containment_a_in_b` instead preferred `Pantala`
- `containment_b_in_a` still preferred `Campodea`

Interpretation of that divergence:

- symmetric overlap-style metrics and directional containment are not answering exactly the same question;
- even on the same `1 Mb` real-sequence extracts, metric choice can alter the top-ranked neighbor;
- the difference here is not dramatic, but it is real and direction-specific.

## Additional pilot pattern

Within this six-taxon masked `1 Mb` subset:

- `Campodea` and `Pantala` were repeatedly near each other across metrics;
- `Catajapyx` was consistently weaker than `Campodea` as a match to `protura_purged`;
- `Folsomia` and `Daphnia` were intermediate but not top-ranked for `protura_purged`.

These are pilot ranking observations only, not phylogenetic conclusions.

## Methodological implication

This pilot suggests that metric choice matters on real sequence content, not only on toy synthetic cases.

What it does support:

- exact full-set overlap metrics and directional containment can disagree on the best match;
- the disagreement is visible even on small real extracts at `k=15`;
- the exact panel is therefore useful as a methodological comparator before altering the production sketch branch.

What it does not support:

- any strong biological interpretation of `protura_purged` affinity;
- any claim that the mini-FASTA ranking should match whole-genome sketch trees;
- any direct replacement of the production high-k workflow.

## Limitations

- only `1 Mb` per taxon was used;
- all mini-FASTA files were truncated at the cap;
- exact sets were held in memory, so this remains a small-pilot tool;
- only masked input was tested here;
- no trees were built;
- no whole-genome exact comparison was attempted;
- results are not directly comparable to Phase 2 full-sketch trees because:
  - input length is capped;
  - only six taxa were used;
  - exact full-set metrics and sketch metrics are different estimators.

## Conservative take

The pilot is sufficient to justify continued methodological work on metric choice.

It is not sufficient to justify biological interpretation.
