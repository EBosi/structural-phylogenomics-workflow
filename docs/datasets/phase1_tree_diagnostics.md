# Phase 1 Small-k Tree Diagnostics

## Scope

This document closes the current small-k assembly-based diagnostics for the
Phase 1 baseline dataset. It uses only outputs already generated under
`results/` and does not introduce new parameters, new runs or new code.

The correct interpretation is conservative: these trees are diagnostics of
small-k genome-wide composition after preprocessing, organelle filtering and
low-complexity masking. They are not a conclusive phylogeny.

## Outputs Produced

- K-mer spectra: `results/kmers/spectra/`
- K-mer feature matrices: `results/kmers/matrices/`
- Distance matrices: `results/distances/`
- Newick trees: `results/trees/`
- Tree manifest: `results/reports/tree_manifest.tsv`
- Pairwise tree comparisons: `results/reports/tree_comparisons.tsv`

Observed output counts:

| Output class | Count |
| --- | ---: |
| per-sample spectra | 40 |
| k-mer matrices | 8 |
| distance matrices | 16 |
| Newick trees | 32 |
| tree manifest rows | 32 |
| tree comparison rows | 496 |

## Parameters Used

All values came from `config/config.yaml`.

| Parameter | Values |
| --- | --- |
| datasets | `unmasked`, `masked` |
| k values | `5`, `6`, `7`, `8` |
| k-mer mode | canonical |
| k-mer normalization | frequency |
| distance metrics | `cosine`, `jensen_shannon` |
| tree methods | `nj`, `upgma` |

## Main Split Frequencies

Across 32 trees, the main internal split counts were:

| Split | Count |
| --- | ---: |
| `Folsomia candida` + `protura_purged` | 15 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 15 |
| `Daphnia pulex` + `protura_purged` | 13 |
| `Daphnia pulex` + `Anopheles gambiae` | 11 |
| `Folsomia candida` + `Anopheles gambiae` | 5 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 3 |
| `Catajapyx aquilonaris` + `Folsomia candida` | 1 |
| `Anopheles gambiae` + `protura_purged` | 1 |

Cherry counts were also mixed:

| Cherry | Count |
| --- | ---: |
| `Daphnia pulex` + `protura_purged` | 13 |
| `Daphnia pulex` + `Anopheles gambiae` | 9 |
| `Folsomia candida` + `protura_purged` | 8 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 5 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 3 |

The recurring `Folsomia candida` + `protura_purged` split is the most plausible
target-near-neighbor signal in this small dataset, but it is not dominant enough
to support a strong phylogenetic claim.

## Stability and Jaccard

Pairwise tree comparison used split-set Jaccard similarity.

| Summary | Value |
| --- | ---: |
| pairwise comparisons | 496 |
| mean Jaccard | 0.309 |
| minimum Jaccard | 0.0 |
| maximum Jaccard | 1.0 |
| comparisons with Jaccard 0.0 | 244 |
| comparisons with Jaccard 0.333 | 148 |
| comparisons with Jaccard 1.0 | 104 |

Grouped comparison means:

| Group | n | Mean Jaccard |
| --- | ---: | ---: |
| same method | 240 | 0.569 |
| cross method | 256 | 0.065 |
| same metric | 240 | 0.328 |
| cross metric | 256 | 0.292 |
| same dataset | 240 | 0.294 |
| cross dataset | 256 | 0.323 |

The largest instability is method-driven: NJ and UPGMA often return different
split sets.

## NJ Versus UPGMA

NJ and UPGMA emphasize different patterns.

NJ split counts:

| Split | Count |
| --- | ---: |
| `Folsomia candida` + `protura_purged` | 15 |
| `Daphnia pulex` + `Anopheles gambiae` | 11 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 5 |
| `Catajapyx aquilonaris` + `Folsomia candida` | 1 |

UPGMA split counts:

| Split | Count |
| --- | ---: |
| `Daphnia pulex` + `protura_purged` | 13 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 10 |
| `Folsomia candida` + `Anopheles gambiae` | 5 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 3 |
| `Anopheles gambiae` + `protura_purged` | 1 |

Interpretation: NJ gives the more biologically plausible
`Folsomia candida` + `protura_purged` signal, while UPGMA appears more sensitive
to global compositional/distance-scale effects. UPGMA should be treated as
diagnostic, not primary evidence.

## Masked Versus Unmasked

Masked split counts:

| Split | Count |
| --- | ---: |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 9 |
| `Folsomia candida` + `protura_purged` | 8 |
| `Daphnia pulex` + `protura_purged` | 7 |
| `Daphnia pulex` + `Anopheles gambiae` | 5 |
| `Folsomia candida` + `Anopheles gambiae` | 2 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 1 |

Unmasked split counts:

| Split | Count |
| --- | ---: |
| `Folsomia candida` + `protura_purged` | 7 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 6 |
| `Daphnia pulex` + `protura_purged` | 6 |
| `Daphnia pulex` + `Anopheles gambiae` | 6 |
| `Folsomia candida` + `Anopheles gambiae` | 3 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 2 |
| `Catajapyx aquilonaris` + `Folsomia candida` | 1 |
| `Anopheles gambiae` + `protura_purged` | 1 |

Masking changes some split frequencies but does not resolve the instability.
The low-complexity mask is useful diagnostically but does not produce a stable
small-k phylogeny by itself.

## Cosine Versus Jensen-Shannon

Cosine split counts:

| Split | Count |
| --- | ---: |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 10 |
| `Folsomia candida` + `protura_purged` | 7 |
| `Daphnia pulex` + `protura_purged` | 7 |
| `Daphnia pulex` + `Anopheles gambiae` | 3 |
| `Folsomia candida` + `Anopheles gambiae` | 3 |
| `Catajapyx aquilonaris` + `Folsomia candida` | 1 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 1 |

Jensen-Shannon split counts:

| Split | Count |
| --- | ---: |
| `Folsomia candida` + `protura_purged` | 8 |
| `Daphnia pulex` + `Anopheles gambiae` | 8 |
| `Daphnia pulex` + `protura_purged` | 6 |
| `Catajapyx aquilonaris` + `Anopheles gambiae` | 5 |
| `Folsomia candida` + `Anopheles gambiae` | 2 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 2 |
| `Anopheles gambiae` + `protura_purged` | 1 |

Both metrics show instability. Jensen-Shannon slightly increases support for
`Folsomia candida` + `protura_purged` and `Daphnia pulex` + `Anopheles gambiae`,
but neither metric yields a robust topology across all k values and methods.

## Nearest-Neighbor Patterns

Nearest-neighbor summaries from distance matrices show that `Daphnia pulex`
is frequently closest to `protura_purged`, especially for k 5-7, despite being
intended as a distant control. This is a warning sign for compositional signal.

Recurring patterns:

- `protura_purged` nearest neighbor is often `Daphnia pulex`, especially across
  k 5-7 and both metrics.
- `Folsomia candida` is often nearest to `protura_purged` at lower k values, but
  shifts toward `Daphnia pulex` or `Catajapyx aquilonaris` at k 8 in some
  conditions.
- `Anopheles gambiae` is often nearest to `Daphnia pulex`, not consistently a
  clean distant outgroup.
- `Catajapyx aquilonaris` alternates among `protura_purged`, `Daphnia pulex` and
  `Anopheles gambiae`, consistent with an unstable signal and possibly affected
  by assembly fragmentation.

This nearest-neighbor behavior is inconsistent with using these small-k trees as
standalone phylogenetic evidence.

## Technical-Biological Interpretation

The current small-k results are best interpreted as a compositional diagnostic:

- k 5-8 frequency spectra are likely dominated by base composition, low-order
  sequence composition, genome-size effects, repeat/low-complexity content and
  assembly properties.
- `protura_purged` has the largest genome in the dataset and the highest masked
  fraction, which can affect frequency-spectrum distances.
- `GCA_000934665.2` is highly fragmented and may carry assembly-quality effects.
- Organelle filtering removed very little sequence overall, but the organelle
  reference set was partial; this can still affect interpretation of any
  organelle-related signal.
- A formal reference tree is required before interpreting concordance or
  discordance biologically.

Conservative conclusion: the baseline remains usable with caution as a
diagnostic dataset, but the small-k assembly trees are too unstable and too
composition-sensitive to support strong phylogenetic claims.
