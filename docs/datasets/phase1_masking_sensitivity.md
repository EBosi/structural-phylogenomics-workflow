# Phase 1 Masking Sensitivity

## Scope

This document evaluates whether the current k-mer/sketch signal is sensitive to
the existing low-complexity/repeat masking. It uses already-produced small-k
outputs and one additional approved micro-target:

- `results/sketch/distances/unmasked/k15/minhash_jaccard.tsv`
- `results/sketch/trees/unmasked/k15/nj.nwk`

No workflow code was modified. No dependency was added. `config/config.yaml` was
not changed. No RepeatMasker, RepeatModeler, EDTA, FASTQ or simulated-read work
was implemented.

## Commands

The first dry-run attempted to use `config/phase1_assembly_baseline.yaml` plus a
temporary sketch-only override. That exposed a real configuration hazard:
`SAMPLE_IDS` was empty because the temporary override replaced top-level
metadata during Snakemake config composition. That command was not used for the
real run.

The corrected dry-run and real run used a single temporary override containing
only the necessary metadata and sketch keys:

```yaml
metadata:
  local_genomes_file: "resources/phase1_assembly_baseline/local_genomes.tsv"
  local_genomes_dir: ""

sketch:
  datasets: ["unmasked"]
  k_values: [15]
  methods: ["nj"]
```

Dry-run:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/sketch/trees/unmasked/k15/nj.nwk \
  results/sketch/distances/unmasked/k15/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile /tmp/<temporary_unmasked_k15_override>.yaml \
  --rerun-triggers mtime
```

Real run:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 \
  results/sketch/trees/unmasked/k15/nj.nwk \
  results/sketch/distances/unmasked/k15/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile /tmp/<temporary_unmasked_k15_override>.yaml \
  --rerun-triggers mtime
```

The valid dry-run planned exactly:

| Rule | Jobs |
| --- | ---: |
| `sketch_signature_per_sample` | 5 |
| `sketch_distance_matrix` | 1 |
| `infer_sketch_tree` | 1 |

The real run completed the same 7 jobs.

## Small-k Masked Versus Unmasked

Inputs used:

- `results/reports/tree_comparisons.tsv`
- `results/distances/masked/`
- `results/distances/unmasked/`
- `results/trees/masked/`
- `results/trees/unmasked/`

Small-k parameters:

| Parameter | Values |
| --- | --- |
| datasets | `masked`, `unmasked` |
| k | `5`, `6`, `7`, `8` |
| metrics | `cosine`, `jensen_shannon` |
| methods | `nj`, `upgma` |

Paired tree comparisons for same k, metric and method:

| Summary | Value |
| --- | ---: |
| paired masked/unmasked comparisons | 16 |
| mean split Jaccard | 0.750000 |
| minimum split Jaccard | 0.333333 |
| maximum split Jaccard | 1.000000 |
| identical split sets | 10 |
| zero-overlap split sets | 0 |

By k:

| k | n | Mean Jaccard | Min | Max |
| --- | ---: | ---: | ---: | ---: |
| 5 | 4 | 1.000000 | 1.000000 | 1.000000 |
| 6 | 4 | 1.000000 | 1.000000 | 1.000000 |
| 7 | 4 | 0.500000 | 0.333333 | 1.000000 |
| 8 | 4 | 0.500000 | 0.333333 | 1.000000 |

By metric and method:

| Group | n | Mean Jaccard | Min | Max |
| --- | ---: | ---: | ---: | ---: |
| cosine | 8 | 0.666667 | 0.333333 | 1.000000 |
| Jensen-Shannon | 8 | 0.833333 | 0.333333 | 1.000000 |
| NJ | 8 | 0.833333 | 0.333333 | 1.000000 |
| UPGMA | 8 | 0.666667 | 0.333333 | 1.000000 |

Interpretation: masking changes some small-k trees, especially at k7-k8, but it
is not the dominant source of instability. The largest instability previously
observed remains method/k/metric dependent rather than a clean masked-versus-
unmasked effect.

## Small-k Split Frequencies

Masked split/cherry counts across 16 small-k trees:

| Split | Count |
| --- | ---: |
| `Anopheles gambiae` + `Catajapyx aquilonaris` | 11 |
| `Folsomia candida` + `protura_purged` | 8 |
| `Daphnia pulex` + `protura_purged` | 7 |
| `Anopheles gambiae` + `Daphnia pulex` | 5 |
| `Anopheles gambiae` + `Folsomia candida` | 2 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 1 |

Unmasked split/cherry counts across 16 small-k trees:

| Split | Count |
| --- | ---: |
| `Folsomia candida` + `protura_purged` | 7 |
| `Anopheles gambiae` + `Catajapyx aquilonaris` | 7 |
| `Daphnia pulex` + `protura_purged` | 6 |
| `Anopheles gambiae` + `Daphnia pulex` | 6 |
| `Anopheles gambiae` + `Folsomia candida` | 3 |
| `Catajapyx aquilonaris` + `Daphnia pulex` | 2 |
| `Catajapyx aquilonaris` + `Folsomia candida` | 1 |
| `Anopheles gambiae` + `protura_purged` | 1 |

Masking reduces the frequency of the suspicious
`Anopheles gambiae` + `Catajapyx aquilonaris` split from 11 to 7, but does not
eliminate instability or produce a single stable topology.

## Small-k Nearest Neighbors

Across the 8 small-k distance matrices per dataset:

| Dataset | `protura_purged` nearest neighbor |
| --- | --- |
| masked | `Daphnia pulex` in 8/8 matrices |
| unmasked | `Daphnia pulex` in 8/8 matrices |

Nearest-neighbor counts for all samples:

| Dataset | Pattern | Count |
| --- | --- | ---: |
| masked | `protura_purged` -> `Daphnia pulex` | 8 |
| masked | `Folsomia candida` -> `protura_purged` | 6 |
| masked | `Anopheles gambiae` -> `Daphnia pulex` | 8 |
| masked | `Catajapyx aquilonaris` -> `Daphnia pulex` | 5 |
| unmasked | `protura_purged` -> `Daphnia pulex` | 8 |
| unmasked | `Folsomia candida` -> `protura_purged` | 5 |
| unmasked | `Anopheles gambiae` -> `Daphnia pulex` | 8 |
| unmasked | `Catajapyx aquilonaris` -> `protura_purged` | 5 |

Masking does not resolve the small-k nearest-neighbor artifact where
`Daphnia pulex` is consistently nearest to `protura_purged`.

## Small-k Distance Shifts

Pairwise small-k distances were compared as:

`unmasked distance - masked distance`

| Summary | Value |
| --- | ---: |
| pairwise comparisons | 80 |
| mean delta | 0.013517 |
| minimum delta | -0.000461 |
| maximum delta | 0.072446 |
| unmasked lower | 1 |
| unmasked higher | 79 |

By metric:

| Metric | n | Mean delta | Unmasked lower | Unmasked higher |
| --- | ---: | ---: | ---: | ---: |
| cosine | 40 | 0.019486 | 1 | 39 |
| Jensen-Shannon | 40 | 0.007548 | 0 | 40 |

By k:

| k | n | Mean delta |
| --- | ---: | ---: |
| 5 | 20 | 0.004283 |
| 6 | 20 | 0.008031 |
| 7 | 20 | 0.014854 |
| 8 | 20 | 0.026900 |

In the small-k frequency-spectrum analysis, unmasked distances are generally
higher, not lower. This does not support a simple interpretation that unmasked
small-k repeats are increasing shared overlap across samples. The effect is
metric- and k-dependent and likely tied to composition/frequency redistribution.

## High-k k15 Masked Versus Unmasked

High-k k15 MinHash/Jaccard outputs compared:

- `results/sketch/distances/masked/k15/minhash_jaccard.tsv`
- `results/sketch/distances/unmasked/k15/minhash_jaccard.tsv`
- `results/sketch/trees/masked/k15/nj.nwk`
- `results/sketch/trees/unmasked/k15/nj.nwk`

All unmasked k15 signatures contain 1000 effective hashes.

Distance summaries:

| Dataset | Min | Max | Mean |
| --- | ---: | ---: | ---: |
| masked | 0.742138 | 0.889506 | 0.817997 |
| unmasked | 0.738170 | 0.881432 | 0.810623 |

High-k k15 nearest neighbors:

| Dataset | Sample | Nearest neighbor | Distance |
| --- | --- | --- | ---: |
| masked | `protura_purged` | `GCA_000934665.2` | 0.742138 |
| unmasked | `protura_purged` | `GCA_000934665.2` | 0.738170 |
| masked | `GCA_000934665.2` | `protura_purged` | 0.742138 |
| unmasked | `GCA_000934665.2` | `protura_purged` | 0.738170 |
| masked | `GCA_021134715.1` | `protura_purged` | 0.853868 |
| unmasked | `GCA_021134715.1` | `protura_purged` | 0.842593 |
| masked | `GCF_943734735.2` | `protura_purged` | 0.774510 |
| unmasked | `GCF_943734735.2` | `protura_purged` | 0.760843 |
| masked | `GCA_052721305.1` | `protura_purged` | 0.783455 |
| unmasked | `GCA_052721305.1` | `protura_purged` | 0.778999 |

Pairwise high-k k15 distance shifts:

| Pair | Masked | Unmasked | Delta |
| --- | ---: | ---: | ---: |
| `protura_purged` - `GCA_000934665.2` | 0.742138 | 0.738170 | -0.003968 |
| `protura_purged` - `GCA_021134715.1` | 0.853868 | 0.842593 | -0.011276 |
| `protura_purged` - `GCF_943734735.2` | 0.774510 | 0.760843 | -0.013667 |
| `protura_purged` - `GCA_052721305.1` | 0.783455 | 0.778999 | -0.004456 |
| `GCA_000934665.2` - `GCA_021134715.1` | 0.889506 | 0.881432 | -0.008074 |
| `GCA_000934665.2` - `GCF_943734735.2` | 0.808105 | 0.798799 | -0.009306 |
| `GCA_000934665.2` - `GCA_052721305.1` | 0.791541 | 0.787144 | -0.004397 |
| `GCA_021134715.1` - `GCF_943734735.2` | 0.871968 | 0.867497 | -0.004471 |
| `GCA_021134715.1` - `GCA_052721305.1` | 0.861048 | 0.851235 | -0.009813 |
| `GCF_943734735.2` - `GCA_052721305.1` | 0.803828 | 0.799520 | -0.004308 |

High-k k15 summary:

| Summary | Value |
| --- | ---: |
| mean delta, unmasked - masked | -0.007374 |
| unmasked lower distances | 10/10 |
| unmasked higher distances | 0/10 |

In contrast to small-k, high-k k15 unmasked distances are consistently lower.
This is compatible with additional shared high-k overlap when masked sequence is
left in, including a possible repeat contribution. It is not by itself proof
that repeats are the biological driver.

## High-k k15 Topology

Masked k15 NJ:

```text
(GCA_052721305.1:0.39244,(GCF_943734735.2:0.39981,(protura_purged:0.35818,GCA_000934665.2:0.38396)Inner1:0.02043)Inner2:0.00757,GCA_021134715.1:0.46861)Inner3:0.00000;
```

Unmasked k15 NJ:

```text
(GCF_943734735.2:0.39539,(protura_purged:0.35493,GCA_000934665.2:0.38324)Inner1:0.01534,(GCA_052721305.1:0.38889,GCA_021134715.1:0.46235)Inner2:0.01250)Inner3:0.00000;
```

Unrooted split-set comparison:

| Comparison | Shared splits | Union splits | Jaccard |
| --- | ---: | ---: | ---: |
| high-k k15 masked vs unmasked NJ | 2 | 2 | 1.000000 |

The Newick display differs, but the unrooted split set is the same:

- `protura_purged` + `GCA_000934665.2`
- `GCA_021134715.1` + `GCA_052721305.1`

Therefore, current masking affects high-k k15 distances but does not change the
unrooted NJ topology for this five-assembly panel.

## Answers

1. How much do masked and unmasked differ in small-k?

They differ modestly. Paired masked/unmasked small-k trees have mean split
Jaccard 0.75, with 10/16 identical split sets and no zero-overlap comparison.
Differences increase at k7-k8, but masking is not the main instability source.

2. Do high-k k15 masked and unmasked give the same nearest neighbor for
`protura_purged`?

Yes. Both give `GCA_000934665.2` (`Catajapyx aquilonaris`) as nearest neighbor
of `protura_purged`.

3. Are unmasked distances lower, suggesting greater overlap due to repeats?

For small-k, no: unmasked distances are higher in 79/80 pairwise comparisons.
For high-k k15, yes: unmasked distances are lower in 10/10 pairwise comparisons.
The high-k result is compatible with additional shared sequence overlap when
masked sequence is retained, but it is not sufficient to attribute the effect
specifically to repeats without a better repeat model and reference comparison.

4. Does masking change topology?

Small-k: sometimes, especially at k7-k8, but not systematically enough to
stabilize the trees. High-k k15: no unrooted topology change; masked and
unmasked NJ have identical split sets.

5. Is there enough evidence to justify a serious repeat backend?

There is enough evidence that masking matters technically, especially for
high-k distance scale. There is not enough evidence that implementing a serious
repeat backend should precede the provisional reference tree. The current
diagnostics cannot yet distinguish repeat-driven signal from broader assembly,
composition and taxon-distance effects.

## Recommendation

Recommendation: B) repeat masking serious useful but after reference tree.

Rationale: the high-k k15 unmasked run shows a consistent distance shift
compatible with repeat/shared-sequence sensitivity, but the topology and
nearest-neighbor result are unchanged. A provisional reference tree is needed
first to decide whether improved repeat masking would correct a phylogenetic
discordance or only perturb already-diagnostic distances.
