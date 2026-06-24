# Phase 1 Reference Hypotheses

## Scope

This document defines three provisional reference hypotheses for the Phase 1
five-assembly panel and records a minimal diagnostic comparison against the
already-produced small-k and high-k NJ trees.

These are reference hypotheses, not a literature-backed reference phylogeny and
not a definitive biological interpretation. They are used only to measure
topological concordance of existing assembly-based outputs.

No workflow code was modified. No new metric was implemented. No FASTQ,
simulated reads, RepeatMasker, RepeatModeler, EDTA or `full_analysis` run was
introduced.

## Taxa

| Sample ID | Label used here |
| --- | --- |
| `protura_purged` | target assembly |
| `GCA_052721305.1` | `Folsomia candida` |
| `GCA_000934665.2` | `Catajapyx aquilonaris` |
| `GCA_021134715.1` | `Daphnia pulex` |
| `GCF_943734735.2` | `Anopheles gambiae` |

## Hypothesis A: `protura_purged` + Folsomia

File:

`docs/datasets/phase1_reference_hypotheses/hypothesis_A_protura_folsomia.nwk`

Newick:

```text
((protura_purged,GCA_052721305.1),GCA_000934665.2,GCA_021134715.1,GCF_943734735.2);
```

What it tests:

This tests whether the assembly-based trees contain the focal split
`protura_purged` + `GCA_052721305.1`.

Assumptions:

It assumes the relevant diagnostic signal would group the target assembly with
the Folsomia assembly before grouping either with the other three samples.

Limits:

This multifurcating topology does not specify relationships among Catajapyx,
Daphnia and Anopheles. It tests one focal split only.

Why provisional:

It is derived from the current diagnostic context, not from a curated reference
tree or literature review. It is therefore only a concordance hypothesis.

Trees expected to be closer if the signal is coherent:

Small-k NJ trees that repeatedly showed `Folsomia candida` + `protura_purged`
should be closest to this hypothesis.

## Hypothesis B: `protura_purged` + Catajapyx

File:

`docs/datasets/phase1_reference_hypotheses/hypothesis_B_protura_catajapyx.nwk`

Newick:

```text
((protura_purged,GCA_000934665.2),GCA_052721305.1,GCA_021134715.1,GCF_943734735.2);
```

What it tests:

This tests whether the assembly-based trees contain the focal split
`protura_purged` + `GCA_000934665.2`.

Assumptions:

It assumes the relevant diagnostic signal would group the target assembly with
the Catajapyx assembly before grouping either with the other three samples.

Limits:

Catajapyx is represented by a fragmented assembly in this baseline. Concordance
with this hypothesis could reflect phylogenetic signal, assembly effects,
repeat/low-complexity effects or other technical structure. This comparison
does not distinguish those causes.

Why provisional:

It is included because the high-k k15 sketch trees currently show this focal
split. That observation is not sufficient to define a biological reference
phylogeny.

Trees expected to be closer if the signal is coherent:

The high-k k15 NJ trees, especially `results/sketch/trees/masked/k15/nj.nwk`
and `results/sketch/trees/unmasked/k15/nj.nwk`, should be closest to this
hypothesis.

## Hypothesis C: Folsomia + Catajapyx

File:

`docs/datasets/phase1_reference_hypotheses/hypothesis_C_folsomia_catajapyx.nwk`

Newick:

```text
((GCA_052721305.1,GCA_000934665.2),protura_purged,GCA_021134715.1,GCF_943734735.2);
```

What it tests:

This tests whether the assembly-based trees group the two non-target
Folsomia/Catajapyx comparator assemblies together instead of grouping either
with `protura_purged`.

Assumptions:

It assumes the relevant diagnostic signal could separate `protura_purged` from
the two comparator assemblies and recover a Folsomia + Catajapyx focal split.

Limits:

This is a control hypothesis for concordance testing. It does not assert that
Folsomia and Catajapyx form a validated biological clade in this panel.

Why provisional:

It is included to test an alternative focal split seen in saturated high-k k21
output. The k21 distance matrix was previously diagnosed as nearly saturated,
so any apparent support from k21 must be declassed.

Trees expected to be closer if the signal is coherent:

Only trees containing the focal split `GCA_052721305.1` + `GCA_000934665.2`
should be closer to this hypothesis. The current k21 tree is not reliable
evidence because it is marked saturated.

## Comparison Method

Output:

`results/reports/phase1_reference_comparison.tsv`

This output is generated and not intended to be tracked.

Compared trees:

| Source | Trees |
| --- | ---: |
| small-k NJ, masked | 8 |
| small-k NJ, unmasked | 8 |
| high-k sketch NJ, k15 masked | 1 |
| high-k sketch NJ, k15 unmasked | 1 |
| high-k sketch NJ, k21 masked | 1 |

Total tree files compared: 19.

Total hypothesis/tree comparisons: 57.

For each comparison, the TSV reports:

| Column | Meaning |
| --- | --- |
| `hypothesis_id` | `A`, `B` or `C` |
| `hypothesis_label` | human-readable focal split |
| `hypothesis_path` | path to the hypothesis Newick file |
| `tree_path` | path to the inferred tree |
| `source_class` | `small_k` or `high_k_sketch` |
| `dataset` | `masked` or `unmasked` |
| `k` | k-mer size |
| `metric` | small-k distance metric or `minhash_jaccard` |
| `method` | tree method, here `nj` |
| `saturated` | `true` for high-k k21, otherwise `false` |
| `shared_splits` | shared internal split count |
| `union_splits` | union internal split count |
| `split_jaccard` | split-set Jaccard |
| `focal_split_present` | whether the tree contains the hypothesis focal split |

Because the hypotheses are multifurcating and encode mainly one focal split,
`focal_split_present` is the primary diagnostic field. `split_jaccard` is a
secondary summary of split-set overlap.

## Command Used

The comparison was generated with an inline Python script that imports
`tree_split_set` and `tree_terminal_names` directly from
`workflow/scripts/phylo_utils.py`.

The script compares:

- `results/trees/masked/*/*/nj.nwk`
- `results/trees/unmasked/*/*/nj.nwk`
- `results/sketch/trees/masked/k15/nj.nwk`
- `results/sketch/trees/unmasked/k15/nj.nwk`
- `results/sketch/trees/masked/k21/nj.nwk`

The high-k masked k21 tree is marked:

```text
saturated=true
```

because the previous high-k sketch report found the k21 distances to be nearly
all `1.000000`.

## Results Summary

Small-k NJ, all 16 trees:

| Hypothesis | Focal split present | Mean split Jaccard |
| --- | ---: | ---: |
| A: `protura_purged` + Folsomia | 15/16 | 0.468750 |
| B: `protura_purged` + Catajapyx | 0/16 | 0.000000 |
| C: Folsomia + Catajapyx | 1/16 | 0.031250 |

Small-k NJ, masked only:

| Hypothesis | Focal split present | Mean split Jaccard |
| --- | ---: | ---: |
| A: `protura_purged` + Folsomia | 8/8 | 0.500000 |
| B: `protura_purged` + Catajapyx | 0/8 | 0.000000 |
| C: Folsomia + Catajapyx | 0/8 | 0.000000 |

Small-k NJ, unmasked only:

| Hypothesis | Focal split present | Mean split Jaccard |
| --- | ---: | ---: |
| A: `protura_purged` + Folsomia | 7/8 | 0.437500 |
| B: `protura_purged` + Catajapyx | 0/8 | 0.000000 |
| C: Folsomia + Catajapyx | 1/8 | 0.062500 |

High-k sketch NJ, k15 only:

| Hypothesis | Focal split present | Mean split Jaccard |
| --- | ---: | ---: |
| A: `protura_purged` + Folsomia | 0/2 | 0.000000 |
| B: `protura_purged` + Catajapyx | 2/2 | 0.500000 |
| C: Folsomia + Catajapyx | 0/2 | 0.000000 |

High-k sketch NJ, k21 saturated:

| Hypothesis | Focal split present | Mean split Jaccard | Status |
| --- | ---: | ---: | --- |
| A: `protura_purged` + Folsomia | 0/1 | 0.000000 | declassed |
| B: `protura_purged` + Catajapyx | 0/1 | 0.000000 | declassed |
| C: Folsomia + Catajapyx | 1/1 | 0.500000 | declassed |

## Interpretation

Small-k NJ is most concordant with Hypothesis A in this diagnostic comparison.
This is a concordance statement only. It is not a biological claim that
Hypothesis A is the true phylogeny.

High-k k15 NJ is most concordant with Hypothesis B. This matches the previous
high-k k15 diagnostic result where `protura_purged` grouped with
`GCA_000934665.2`.

High-k k21 is excluded from support ranking because it is saturated. Although it
contains the Hypothesis C focal split, the k21 distance matrix was previously
diagnosed as nearly non-informative for this panel.

The discordance between small-k NJ and high-k k15 should not be interpreted as
biological without a more formal external reference and broader taxon sampling.

## Conservative Next Step

The next step should be a formalized reference-tree comparison layer or a
manually curated provisional reference based on external taxonomy/literature,
kept separate from the current diagnostic hypotheses. FASTQ, simulated reads and
serious repeat masking should remain out of scope until this reference
comparison is clarified.
