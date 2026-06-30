# Phase 2 High-k Sketch Diagnostics

## Scope

This document summarizes the minimal Phase 2 high-k sketch run on the 16-sample expanded assembly panel.

Run scope:

- datasets: `masked`, `unmasked`
- k: `15`, `21`
- sketch distance: `minhash_jaccard`
- tree method: `nj`

Guardrails respected:

- no FASTQ analysis
- no read simulation
- no repeat-backend implementation
- no frequency filtering
- no `k17`
- no `upgma`
- no writes outside `results/phase2_expanded_panel/` and `resources/phase2_expanded_panel/`

## Commands

Dry-run:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/phase2_expanded_panel/sketch/trees/masked/k15/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/unmasked/k15/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/masked/k21/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/unmasked/k21/nj.nwk \
  results/phase2_expanded_panel/sketch/distances/masked/k15/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/unmasked/k15/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/masked/k21/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/unmasked/k21/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

Run:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 \
  results/phase2_expanded_panel/sketch/trees/masked/k15/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/unmasked/k15/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/masked/k21/nj.nwk \
  results/phase2_expanded_panel/sketch/trees/unmasked/k21/nj.nwk \
  results/phase2_expanded_panel/sketch/distances/masked/k15/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/unmasked/k15/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/masked/k21/minhash_jaccard.tsv \
  results/phase2_expanded_panel/sketch/distances/unmasked/k21/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

## Outputs

- sketch signatures: `64`
- sketch distance matrices: `4`
- NJ trees: `4`
- failed dataset/k combinations: `0`

All `64/64` sketch signatures contain `1000` retained hashes, matching the configured `num_hashes`.

Produced targets:

- `results/phase2_expanded_panel/sketch/distances/masked/k15/minhash_jaccard.tsv`
- `results/phase2_expanded_panel/sketch/distances/unmasked/k15/minhash_jaccard.tsv`
- `results/phase2_expanded_panel/sketch/distances/masked/k21/minhash_jaccard.tsv`
- `results/phase2_expanded_panel/sketch/distances/unmasked/k21/minhash_jaccard.tsv`
- `results/phase2_expanded_panel/sketch/trees/masked/k15/nj.nwk`
- `results/phase2_expanded_panel/sketch/trees/unmasked/k15/nj.nwk`
- `results/phase2_expanded_panel/sketch/trees/masked/k21/nj.nwk`
- `results/phase2_expanded_panel/sketch/trees/unmasked/k21/nj.nwk`

## Saturation Diagnostics

Pairwise distance counts are over `120` unordered pairs for `16` taxa.

### masked k15

- min distance: `0.5663`
- max distance: `0.9230`
- mean distance: `0.7915`
- pairs `>= 0.99`: `0/120`
- pairs `== 1.0`: `0/120`
- classification: `informative`

### unmasked k15

- min distance: `0.5724`
- max distance: `0.9148`
- mean distance: `0.7853`
- pairs `>= 0.99`: `0/120`
- pairs `== 1.0`: `0/120`
- classification: `informative`

### masked k21

- min distance: `0.9960`
- max distance: `1.0000`
- mean distance: `0.9994`
- pairs `>= 0.99`: `120/120`
- pairs `== 1.0`: `46/120`
- classification: `saturated`

### unmasked k21

- min distance: `0.9919`
- max distance: `1.0000`
- mean distance: `0.9984`
- pairs `>= 0.99`: `120/120`
- pairs `== 1.0`: `7/120`
- classification: `saturated`

Bottom line on saturation:

- `k15` is usable as a minimal high-k diagnostic
- `k21` is not reliable for biological interpretation in this panel and should remain a saturation control only

## Nearest-Neighbor Diagnostics

### protura_purged

Primary k15 result:

- `masked k15`: nearest neighbor `GCA_020796165.1` / *Pantala flavescens*
- `unmasked k15`: nearest neighbor `GCA_020796165.1` / *Pantala flavescens*

Additional k15 context:

- `Campodea augens` is also close to `protura_purged`
- in `masked k15`, `Campodea`, `Pantala`, and `Jaera` are tied at the same minimum distance to `protura_purged`
- in `unmasked k15`, `Campodea` is second closest, behind `Pantala`
- `Catajapyx aquilonaris` is not the nearest neighbor in either `k15` dataset
- `Folsomia candida` is clearly farther than both diplurans in `k15`

Saturated control:

- `masked k21`: nearest neighbor `GCA_947179485.1` / *Allacma fusca*
- `unmasked k21`: nearest neighbor `GCA_009757345.1` / *Campodea augens*

These `k21` nearest-neighbor calls are not trustworthy because the matrices are saturated.

### Diplura

- `Catajapyx aquilonaris` nearest neighbor:
  - `masked k15`: `GCA_947179485.1` / *Allacma fusca*
  - `unmasked k15`: `GCA_947179485.1` / *Allacma fusca*
- `Campodea augens` nearest neighbor:
  - `masked k15`: `GCA_020796165.1` / *Pantala flavescens*
  - `unmasked k15`: `GCA_020796165.1` / *Pantala flavescens*

So the two Diplura do not preferentially retrieve each other at `k15`.

### Diplura-vs-Protura patterns

- `protura_purged` is closer to `Campodea augens` than to `Catajapyx aquilonaris` in both `k15` datasets
- `protura_purged` is also closer to `Pantala flavescens` than to either dipluran in both `k15` datasets
- `Catajapyx` and `Campodea` do not form a direct cherry in any of the four trees

### Collembola

- `protura_purged` does not take any collembolan as nearest neighbor in `k15`
- `protura_purged` does not form a direct cherry with a collembolan in any of the four high-k trees
- `Folsomia candida` is not favored over Diplura at `k15`

## Topology Diagnostics

### masked k15

Main local pattern around `protura_purged`:

- `Campodea augens + Pantala flavescens` is a direct cherry
- that cherry joins `Thermobia domestica`
- `protura_purged` joins that mixed clade next

Implication:

- this is not a direct `Protura + Diplura` cherry
- the local neighborhood is mixed and includes one dipluran plus ectognath insects

Other notable cherries:

- `Folsomia candida + Lepidocyrtus curvicollis`
- `Tigriopus californicus + Cloeon dipterum`
- `Daphnia pulex + Podura aquatica`

### unmasked k15

Main local pattern around `protura_purged`:

- `Campodea augens + Pantala flavescens` remains a direct cherry
- that cherry joins `Thermobia domestica`
- `protura_purged` again joins this mixed clade
- `Allacma fusca` then joins next
- `Catajapyx aquilonaris` is outside that local core

Other notable cherries:

- `Orchesella flavescens + Lepidocyrtus curvicollis`
- `Jaera ischiosetosa + Pogonognathellus longicornis`
- `Daphnia pulex + Podura aquatica`

### masked vs unmasked at k15

Shared behavior:

- both datasets place `protura_purged` in a mixed local neighborhood containing `Campodea`, `Pantala`, and `Thermobia`
- both datasets avoid the Phase 1 `protura_purged + Catajapyx` signal
- both datasets avoid the Phase 1 small-k `protura_purged + Folsomia` signal

Differences:

- `masked k15` pulls `Catajapyx` later and keeps a `Folsomia + Lepidocyrtus` cherry
- `unmasked k15` gives a somewhat different collembolan/insect arrangement, but the local `protura` neighborhood is broadly similar

### k21 control

`k21` trees are not biologically interpretable here.

Reasons:

- all pairwise distances are `>= 0.99`
- branch lengths are nearly flat around `0.498-0.500`
- cherries differ between `masked` and `unmasked` despite pervasive saturation

The only safe use of `k21` in Phase 2 is to document that this sketch setting is too saturated for this panel.

## Comparison to Phase 1

Phase 1 high-k `k15` supported:

- `protura_purged + Catajapyx aquilonaris`

Phase 2 high-k `k15` does not preserve that pattern.

Observed instead:

- `protura_purged` is nearest to `Pantala flavescens` in both `k15` matrices
- `Campodea augens` is closer to `protura_purged` than `Catajapyx aquilonaris` in both `k15` matrices
- `Catajapyx` is not the nearest neighbor of `protura_purged`
- the two Diplura do not define a stable dipluran-only local cluster around `protura_purged`

Conservative interpretation:

- adding a second dipluran does not recover a clear Diplura-wide high-k signal
- the original Phase 1 `protura_purged + Catajapyx` result looks taxon-sampling-sensitive and may have been Catajapyx-specific or panel-limited

`k21` remains saturated in Phase 2, consistent with the Phase 1 decision to declass it.

## Comparison to Phase 2 Small-k

Phase 2 small-k result:

- `protura_purged` nearest to `Daphnia pulex` in `16/16` matrices

Phase 2 high-k `k15` changes that pattern.

Observed change:

- `Daphnia pulex` is no longer the nearest neighbor of `protura_purged`
- `k15` places `protura_purged` closer to a mixed neighborhood containing `Pantala`, `Campodea`, and `Thermobia`

This is more interpretable than the small-k `Protura -> Daphnia` attraction, but it is still not a clean phylogenetic signal:

- the local pattern is mixed rather than taxonomically coherent
- no direct `Protura + Diplura` cherry is recovered
- no direct `Protura + Collembola` cherry is recovered

So high-k is preferable to small-k here, but still not decisive.

## Caveats

- no strong biological claims should be made from this run
- assembly size varies strongly across the panel
- `Thermobia domestica` is extremely large
- `Catajapyx aquilonaris`, `Campodea augens`, and `protura_purged` are fragmented
- mitochondrial references are incomplete for several key taxa
- sketch size is limited to `1000` retained hashes
- `k17` was intentionally skipped to save time and can be tested later only if it is needed to bridge `k15` and saturated `k21`
- frequency filtering remains a future idea and was not implemented here

## Conservative Conclusion

Headline high-k result:

- Phase 2 high-k `k15` does not preserve the Phase 1 `protura_purged + Catajapyx` signal.
- Instead, `protura_purged` sits in a mixed local neighborhood led by `Pantala flavescens`, with `Campodea augens` closer than `Catajapyx aquilonaris`.
- `k21` is saturated in both masked and unmasked runs and should remain declassed.

Conservative next-step judgment:

- `k17` is not required immediately.
- The current evidence is sufficient to proceed to a formal Phase 2 signal-comparison step, because:
  - small-k is clearly unstable and compositional
  - high-k `k15` is more informative than small-k
  - high-k `k15` still does not yield a clean taxonomic resolution
  - `k21` is confirmed as saturated

That next comparison should remain technical and reference-based rather than strongly biological.
