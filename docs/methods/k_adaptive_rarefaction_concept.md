# k-adaptive rarefaction concept

## Purpose

This document outlines a possible new direction for the workflow: a multi-k, sampling-aware alignment-free distance framework based on k-mer response profiles.

The motivation is that a single fixed value of `k` is unlikely to capture robust signal across divergent, fragmented, repeat-rich, and unevenly assembled genomes.

Instead of asking:

> Which single k-mer size should we use?

this approach asks:

> How does the k-mer representation of a genome, and the similarity between genomes, behave as k changes?

The central idea is to treat `k` as a scale parameter and to estimate, for each `k`, how much sequence must be sampled to obtain a stable representation of a genome.

## Core idea

For each k-mer size `k`, define a sampling effort:

```text
L = f(k)
```

where `L(k)` is the amount of genomic sequence needed to obtain a sufficiently stable k-mer representation at that value of `k`.

Small values of `k` may require less sequence because the k-mer space is small and saturates quickly.

Larger values of `k` may require more sequence because the k-mer space is larger and local sampling variation becomes more important.

The goal is not to compare whole genomes directly at first. The first goal is to estimate stable genome-level k-mer behavior through standardized replicated sampling.

## Key distinction

This approach separates two sources of variation:

1. **Intra-genome sampling variation**

   How much do k-mer spectra or k-mer sets vary when different regions of the same genome are sampled?

2. **Inter-genome difference**

   Once intra-genome sampling variation is understood, how different are two genomes from each other?

The calibration of `L(k)` should be based on intra-genome stability, before comparing different genomes.

## Step 1 — Intra-genome calibration

For a genome `G`, a k-mer size `k`, and a sampling length `L`:

1. draw multiple replicated samples from `G`;
2. compute a k-mer representation for each replicate;
3. compare replicate representations within the same genome;
4. measure how stable the representation is as `L` increases.

Conceptually:

```text
G_sample_1(k, L) -> spectrum/set/sketch 1
G_sample_2(k, L) -> spectrum/set/sketch 2
...
G_sample_n(k, L) -> spectrum/set/sketch n
```

Then estimate the distribution of intra-genome distances:

```text
D_intra(G, k, L)
```

The purpose is not necessarily to make different genomic regions identical. For larger `k`, non-overlapping regions of the same genome may remain different. The goal is to find when the distribution of intra-genome distances becomes stable as `L` increases.

## Step 2 — Defining L(k)

For each `k`, choose a global sampling length `L(k)` based on intra-genome stability.

Possible criteria:

```text
L(k) = smallest L where increasing L further gives little improvement
```

or:

```text
L(k) = smallest L where intra-genome distance variance reaches a plateau
```

or:

```text
L(k) = smallest L that stabilizes most genomes in the panel
```

At this stage, `L(k)` should not depend on the pair of genomes being compared.

A conservative version is:

```text
L(k) = value that stabilizes 80–90% of genomes in the panel
```

or, for a stricter pilot:

```text
L(k) = maximum required L among the focal genomes
```

## Step 3 — Inter-genome comparison

After `L(k)` has been chosen, compare genomes using the same sampling effort at each `k`.

For two genomes `A` and `B`:

```text
D_inter(A, B, k, L(k))
```

This can then be interpreted relative to intra-genome sampling variation:

```text
D_inter(A, B, k, L(k))
vs
D_intra(A, k, L(k)) and D_intra(B, k, L(k))
```

A useful normalized quantity could be:

```text
D_normalized(A, B, k) =
D_inter(A, B, k, L(k)) - expected_intra_variation(A, B, k)
```

or a z-score-like quantity:

```text
Z(A, B, k) =
(D_inter - mean_intra) / sd_intra
```

This would help distinguish true inter-genome separation from sampling noise.

## Step 4 — Multi-k response profile

For each genome pair, compute a profile across multiple k values:

```text
D(A, B, k=7)
D(A, B, k=9)
D(A, B, k=11)
D(A, B, k=13)
D(A, B, k=15)
D(A, B, k=17)
D(A, B, k=19)
...
```

The shape of this curve is the main signal.

Possible interpretations:

* similarity only at low `k` may indicate compositional similarity;
* persistence of similarity into mid/high `k` may indicate more specific shared sequence content;
* strong unmasked but weak masked signal may indicate repeat-driven similarity;
* unstable curves may indicate insufficient sampling, assembly artifacts, or weak signal.

## Candidate outputs

The method may produce several diagnostic quantities rather than a single distance at first:

```text
mean similarity across k
decay slope across k
area under the similarity curve
sampling stability
masking sensitivity
highest k with detectable signal
```

A possible final distance could be based on:

```text
1 - AUC of the normalized similarity curve
```

with penalties for:

```text
high sampling instability
high masking sensitivity
lack of plateau at larger k
```

## Sampling strategy

Initial sampling should be standardized and reproducible.

Suggested strategy:

* split each genome into non-overlapping windows;
* sample windows without replacement within each replicate;
* allow windows to be reused across different replicates;
* use a fixed random seed;
* preserve source coordinates in a manifest.

Initial window size:

```text
50 kb
```

This is a compromise between local genomic structure and enough granularity for random sampling.

## Candidate k and L grid for calibration

Before defining a final `L(k)`, test a grid such as:

```text
k = 7, 9, 11, 13, 15, 17, 19
L = 0.25 Mb, 0.5 Mb, 1 Mb, 2 Mb, 5 Mb, 10 Mb, 20 Mb
```

The goal of the grid is not to use all combinations permanently. The goal is to estimate which `L` is sufficient for each `k`.

A possible later adaptive schedule could look like:

```text
k=7   -> L=0.25 Mb
k=9   -> L=0.5 Mb
k=11  -> L=1 Mb
k=13  -> L=2 Mb
k=15  -> L=5 Mb
k=17  -> L=10 Mb
k=19  -> L=20 Mb
```

These values are placeholders and must be empirically calibrated.

## Representation choices

Possible k-mer representations:

### 1. Frequency spectra

Pros:

* normalized by sequence length;
* useful for low and mid k;
* captures compositional signal.

Cons:

* high-dimensional and sparse at larger k;
* may be dominated by composition at low k.

### 2. Unique k-mer sets

Pros:

* simple;
* compatible with Jaccard-like comparisons;
* conceptually close to current high-k sketch work.

Cons:

* sensitive to sample length;
* not scalable for large exact comparisons;
* loses abundance information.

### 3. Sketches

Pros:

* scalable;
* closer to eventual production use.

Cons:

* current fixed-size bottom-k sketches are not suitable for true containment;
* sketch design must match the intended metric.

The first prototype should probably use exact sets or frequency spectra on sampled windows, because the goal is methodological calibration rather than production speed.

## Metrics to consider

For intra-genome stability:

```text
Jaccard distance between sampled k-mer sets
cosine distance between frequency spectra
Jensen-Shannon distance between frequency spectra
number of unique k-mers
variance in unique-kmer count
```

For inter-genome comparison:

```text
Jaccard / Mash-like distance
Dice similarity
cosine similarity
frequency-based distance
normalized distance relative to intra-genome variation
```

Containment should be treated separately, because it answers a directional query/reference question rather than a symmetric taxon-distance question.

## Relationship to previous work

Previous diagnostics showed that:

* small-k spectra can be dominated by compositional signal;
* k15 bottom-k sketches produce usable but masking-sensitive signal;
* k21 sketches were largely saturated or uninformative;
* exact mini-FASTA panels showed that metric choice can matter;
* fixed-size bottom-k signatures do not support independent containment-style metrics because all signatures have the same retained hash count.

This motivates a shift from single-k distances to multi-k response profiles.

## First minimal pilot

A minimal pilot should not attempt a full phylogeny.

Suggested pilot:

* focal genome: `protura_purged`;
* comparators:

  * `Campodea`
  * `Catajapyx`
  * `Folsomia`
  * `Pantala`
  * `Daphnia`
* input: masked FASTA first;
* k values:

  * `7, 9, 11, 13, 15, 17, 19`
* L values:

  * `0.25, 0.5, 1, 2, 5, 10, 20 Mb`
* replicates:

  * start with `10–20`;
* no trees;
* no Snakemake integration initially;
* output tables and diagnostic plots only.

Primary questions:

1. For each `k`, how large must `L` be before intra-genome k-mer behavior stabilizes?
2. Does `L(k)` increase with `k` as expected?
3. Do inter-genome differences exceed intra-genome sampling variation?
4. Which genome pairs maintain similarity across increasing `k`?
5. Which signals are masking-sensitive?

## Expected outcomes

Useful outcome:

```text
Stable intra-genome behavior emerges at increasing L as k increases.
Inter-genome profiles differ more than intra-genome sampling variation.
Some genome pairs show persistent mid-k similarity.
```

Problematic but informative outcome:

```text
Intra-genome variation remains high even at large L.
Inter-genome differences are not clearly separable from sampling variation.
```

This would suggest that the available assemblies or k-mer representation are not sufficient for robust distance estimation at the tested scale.

## Current decision

This concept should be developed as a new branch of the project.

The previous k-mer distance diagnostics should remain available as methodological background, but further work should focus on:

```text
multi-k response
k-adaptive sampling length
intra-genome calibration
inter-genome normalization
```

rather than adding more formulas to the current fixed-size bottom-k sketch workflow.
