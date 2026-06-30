# Phase 2 Small-k NJ Diagnostics

## Scope

This document summarizes the Phase 2 assembly-only small-k NJ run on the 16-sample expanded panel under `results/phase2_expanded_panel/`.

Run scope:

- datasets: `masked`, `unmasked`
- k: `5`, `6`, `7`, `8`
- metrics: `cosine`, `jensen_shannon`
- method: `nj` only

No FASTQ analysis, high-k sketch, UPGMA, frequency filtering, or new masking backends were used here.

## Commands

Dry-run:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/phase2_expanded_panel/reports/tree_manifest.tsv \
  results/phase2_expanded_panel/reports/tree_comparisons.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

Run:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 \
  results/phase2_expanded_panel/reports/tree_manifest.tsv \
  results/phase2_expanded_panel/reports/tree_comparisons.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

## Output Counts

- k-mer spectra: `128`
- k-mer feature matrices: `8`
- distance matrices: `16`
- NJ trees: `16`
- UPGMA trees: `0`

UPGMA was successfully avoided.

## Dataset Sanity

- Samples in `pre_kmer_summary.tsv`: `16`
- Expected samples present: yes

Known caveats carried from Phase 2 pre-kmer:

- `GCA_964235325.1` / *Thermobia domestica* is extremely large, about `5.13 Gbp`
- `GCA_000934665.2` / *Catajapyx aquilonaris* is highly fragmented
- `GCA_009757345.1` / *Campodea augens* is large and fragmented
- `protura_purged` is fragmented
- `GCA_046254475.1` / *Podura aquatica* has high masking fraction
- organism-specific mitochondrial references are missing for several key taxa, including `protura_purged`, `Catajapyx`, and `Campodea`

## Small-k Topology Diagnostics

### Protura nearest-neighbor from distance matrices

`protura_purged` has the same nearest neighbor in all `16/16` small-k distance matrices:

- nearest neighbor: `GCA_021134715.1` / *Daphnia pulex*
- broad group: `crustacea`

This is the dominant small-k signal in Phase 2.

Secondary ranking is parameter-sensitive:

- `Folsomia candida` is often rank `2` at lower k, especially `k5-k6`
- `Catajapyx aquilonaris` becomes relatively closer in some unmasked higher-k small-k settings, especially `unmasked k7-k8`
- `Campodea augens` is never the nearest neighbor and usually ranks low

### Protura placement in NJ trees

There is no consistent `protura_purged + Collembola` or `protura_purged + Diplura` cherry across the 16 trees.

Observed patterns:

- `15/16` trees place `protura_purged` next to a mixed sibling clade containing some combination of crustacean, insect, and dipluran taxa
- `1/16` tree, `unmasked k8 cosine`, gives a direct `protura_purged + Campodea augens` cherry
- `0/16` trees give a direct `protura_purged + Folsomia candida` cherry
- `0/16` trees give any immediate Collembola-only sister grouping to `protura_purged`

### Group coherence

Across these `16` NJ trees:

- Collembola are not monophyletic in any tree
- Diplura are not monophyletic in any tree
- Crustacea are not monophyletic in any tree

This means the expanded Phase 2 small-k trees are not recovering clean internal coherence for these broad groups.

### Crustacea as outgroups

Crustacea do not behave as clean distant outgroups in small-k Phase 2.

Instead:

- `Daphnia pulex` is the nearest neighbor of `protura_purged` in all distance matrices
- crustacean taxa are repeatedly pulled into the local neighborhood of `protura_purged`

That behavior is more consistent with compositional attraction or other non-phylogenetic small-k effects than with a trustworthy deep outgroup pattern.

## Parameter Sensitivity

### Masked vs unmasked

- `masked` trees are more internally stable than `unmasked` trees
- `unmasked cosine` is the least stable combination, especially at higher small-k
- `masked` runs still show the same broad problem: attraction of `protura_purged` toward mixed clades involving crustacean and non-collembolan taxa

### k5, k6, k7, k8

- lower-k and mid-k runs do not recover a robust Phase 1-like `Protura + Folsomia` tree signal
- `unmasked k8 cosine` is the most deviant setting and produces the only direct `protura_purged + Campodea augens` cherry
- even where topology changes, the distance-level nearest neighbor of `protura_purged` remains `Daphnia pulex`

### Cosine vs Jensen-Shannon

Using split-set Jaccard from `tree_comparisons.tsv`:

- `masked cosine` is highly stable across adjacent k values, with Jaccard about `0.86-1.00`
- `masked jensen_shannon` is also stable, with Jaccard about `0.86-1.00`
- `unmasked jensen_shannon` is moderately stable and sometimes identical across non-adjacent k
- `unmasked cosine` is clearly unstable, dropping to Jaccard `0.44` between `k6-k7` and `0.08` between `k7-k8`

So the strongest instability is not simply "all small-k", but specifically the unmasked cosine branch at higher small-k.

## Phase 1 Comparison

Phase 1 small-k NJ mostly supported hypothesis A:

- `protura_purged + Folsomia`

Phase 2 small-k does not preserve that signal as a tree-level result.

What remains from the Phase 1 pattern:

- `Folsomia candida` is often among the relatively closer taxa to `protura_purged` at lower k

What changes in Phase 2:

- the dominant nearest-neighbor signal is `Daphnia pulex`, not a collembolan
- tree topology does not support a stable Ellipura-like `Protura + Collembola` pattern
- only one parameter set yields a direct dipluran cherry, `protura_purged + Campodea`, and that is an isolated unstable case

Conservative reading:

- expanded taxon sampling makes the Phase 1 small-k signal look unstable rather than reinforced
- Phase 2 small-k does not provide robust support for either a Protura+Collembola or Protura+Diplura interpretation

## Conservative Conclusion

- No strong biological or phylogenetic claims should be made from these Phase 2 small-k NJ trees.
- Small-k signal remains plausibly compositional.
- Expanded sampling did not stabilize the Phase 1 small-k pattern.
- The strongest small-k attraction, `protura_purged -> Daphnia pulex`, is difficult to interpret biologically and should be treated as suspect.
- High-k Phase 2 sketch analysis is still needed, and is now the appropriate next computational step.

## Bottom Line

Headline small-k result:

- Phase 2 small-k NJ is dominated by unstable non-reference-like signal, with `protura_purged` consistently nearest to `Daphnia pulex` in distance space and no stable recovery of a Protura+Collembola tree pattern.
