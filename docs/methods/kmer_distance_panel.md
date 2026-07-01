# k-mer distance panel

## Purpose

This note defines a small validation panel for set-based k-mer distances before further biological interpretation of the high-k sketch branch.

The current production sketch workflow uses:

- unique k-mer presence/absence;
- canonical k-mers;
- bottom-k / smallest-hash sketches;
- Jaccard distance on retained sketch hash sets.

That is a reasonable fast approximation to a symmetric overlap metric, but it does not by itself address:

- genome-size imbalance;
- subset-like relationships;
- partial or fragmented assemblies;
- future read-like or low-coverage data.

The validation panel therefore uses exact unique-kmer sets on synthetic cases to compare several metrics in a controlled setting.

## Metrics in the panel

For k-mer sets `A` and `B`:

- `n_a = |A|`
- `n_b = |B|`
- `intersection = |A ∩ B|`
- `union = |A ∪ B|`

Derived metrics:

- `jaccard = intersection / union`
- `jaccard_distance = 1 - jaccard`
- `dice = 2 * intersection / (n_a + n_b)`
- `dice_distance = 1 - dice`
- `overlap_coefficient = intersection / min(n_a, n_b)`
- `overlap_distance = 1 - overlap_coefficient`
- `containment_a_in_b = intersection / n_a`
- `containment_b_in_a = intersection / n_b`
- `max_containment = intersection / max(n_a, n_b)`
- `binary_cosine = intersection / sqrt(n_a * n_b)`
- `binary_cosine_distance = 1 - binary_cosine`
- `mash_distance = -(1 / k) * ln(2J / (1 + J))` for `J > 0`; `inf` if `J = 0`

## Symmetry and likely use

Symmetric metrics:

- Jaccard
- Dice
- overlap coefficient
- max-containment
- binary cosine
- Mash distance

These are easier to use in tree-building pipelines because they define a single pairwise quantity for `A,B`.

Directional metrics:

- containment `A in B`
- containment `B in A`

These are often better for placement or query-vs-reference questions, because they can distinguish:

- “how much of the small/partial query is found in the reference?”
- from
- “how much of the full reference is recovered by the query?”

That distinction is exactly what symmetric Jaccard suppresses.

## What each metric measures

### Jaccard

Measures shared content relative to total union.

Useful for:

- complete sets of roughly comparable size and completeness;
- exact set-overlap baselines;
- fast symmetric screening.

Weakness:

- penalizes large private tails;
- penalizes genome-size imbalance;
- treats subset-like relationships as weak similarity if the larger set has many extra k-mers.

### Dice

Measures shared content relative to total size, but is less punitive than Jaccard.

Useful for:

- a softer symmetric overlap metric;
- cases where Jaccard seems too harsh.

Weakness:

- still symmetric;
- still declines when one set has a large private tail;
- still does not represent directional containment.

### Overlap coefficient

Measures shared content relative to the smaller set.

Useful for:

- asking whether the smaller set is fully represented in the larger one;
- subset-like comparisons in a symmetric wrapper.

Weakness:

- can report `1.0` even when one genome has a very large private tail;
- may overstate “closeness” if the real question is global similarity rather than recovery of the smaller set.

### Containment

Measures directional recovery:

- `A in B` asks how much of `A` is represented in `B`;
- `B in A` asks the opposite.

Useful for:

- partial assemblies;
- future read-level data;
- target-vs-reference placement;
- strong genome-size imbalance.

Weakness:

- directional, so less natural for standard tree construction;
- does not directly answer “global similarity.”

### Max-containment

Measures overlap relative to the larger set.

Useful for:

- a stricter symmetric quantity when the larger set should dominate interpretation.

Weakness:

- very punitive to partial or query-like data;
- rarely the best summary for subset-like problems.

### Binary cosine

Measures angular similarity of binary membership vectors.

Useful for:

- a symmetric measure less harsh than Jaccard in some imbalanced cases;
- comparisons where some normalization for set size is desirable.

Weakness:

- still symmetric;
- still does not solve the directional placement problem.

### Mash distance

Mash distance is a Jaccard-derived transform, not a new solution to containment mismatch.

Useful for:

- complete-ish genome comparisons under a divergence-style interpretation;
- a familiar alignment-free quantity when Jaccard behavior is already acceptable.

Weakness:

- inherits the Jaccard dependence on global set overlap;
- does not solve subset/containment asymmetry;
- does not become read-aware simply by using sketches.

## Why Jaccard penalizes private content

Jaccard uses:

- numerator: shared content
- denominator: shared content plus all private content

So when one set has many extra k-mers:

- shared core may stay constant;
- union grows;
- Jaccard drops.

This makes Jaccard sensitive to:

- genome-size imbalance;
- assembly completeness imbalance;
- fragmented target vs larger reference;
- reads or low-coverage data with many missing true k-mers.

That behavior is mathematically correct for Jaccard. The question is whether it matches the biological question.

## Why containment may fit partial/read-like problems better

Containment asks:

- “how much of the query is found in the reference?”

That is often closer to the intended question for:

- fragmented assemblies;
- incomplete target genomes;
- low-coverage or read-like data;
- query-vs-reference workflows.

Containment does not pretend that “query missing most of the reference” means the query is biologically distant. It separates:

- incomplete recovery
- from
- global lack of shared content.

## Why exact Jaccard still matters

Exact Jaccard remains useful as a baseline because it lets us distinguish:

- estimator problems from bottom-k sketching;
- from
- conceptual problems inherent to Jaccard itself.

If exact Jaccard already behaves poorly on a toy subset case, then improving the sketch estimator alone will not solve the biological mismatch.

## Zero-size behavior in the validation panel

The test-only validation panel treats undefined denominators explicitly:

- metrics with denominator `0` are returned as `None`
- Mash distance is `None` if Jaccard is undefined
- Mash distance is `inf` if Jaccard is exactly `0`

This is intentionally explicit. It avoids silently pretending that empty-vs-empty or empty-query behavior is biologically meaningful.

## Main methodological implication

Exact Jaccard is a useful baseline.

It is not necessarily the right final metric for:

- partial assemblies;
- fragmented targets vs larger references;
- future read-compatible workflows.

Containment-like quantities are therefore likely to matter if the long-term goal is robust query/reference comparison rather than only symmetric complete-genome trees.
