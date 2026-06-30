# Phase 2 Signal Comparison

## Phase 2 Objective

Phase 2 tested whether the two Phase 1 signals persisted after expanding the assembly-only taxon panel from 5 taxa to 16 taxa.

Phase 1 anchor observations were:

- small-k suggested `protura_purged + Folsomia candida`
- high-k `k15` suggested `protura_purged + Catajapyx aquilonaris`
- high-k `k21` was saturated and declassed

The Phase 2 goal was not to prove a phylogeny. It was to ask whether either Phase 1 signal survived broader taxon sampling, or whether both signals broke down under a more demanding panel.

## Dataset Context

Phase 2 used the expanded 16-assembly panel documented in [phase2_expanded_panel.md](/media/shared1/structural-phylogenomics-workflow/docs/datasets/phase2_expanded_panel.md).

Relevant caveats from existing files:

- `protura_purged` is fragmented
- `Catajapyx aquilonaris` is highly fragmented
- `Campodea augens` is large and fragmented
- `Thermobia domestica` is extremely large, about `5.13 Gbp`
- `Podura aquatica` has a high masking fraction
- mitochondrial references are incomplete for several key taxa
- assembly size and fragmentation vary strongly across the panel

These caveats matter because the current workflow is being evaluated as a signal-diagnostic framework, not yet as a finalized tree-first phylogenetic method.

## Summary Table

| Signal class | Informativeness | `protura_purged` nearest-neighbor behavior | Masked/unmasked sensitivity | Supports Protura+Collembola | Supports Protura+Diplura | Taxonomic coherence | Use for downstream method development |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Phase 2 small-k `k5-k8` | Low for phylogenetic interpretation; compositionally unstable | `Daphnia pulex` in `16/16` distance matrices | Moderate; `unmasked cosine` most unstable, but core `Protura -> Daphnia` effect persists | No robust support | No robust support; one isolated `protura_purged + Campodea` tree only | Poor; Collembola, Diplura, and Crustacea are not coherent in the NJ trees | Useful as a negative-control diagnostic for compositional instability, not as a primary tree signal |
| Phase 2 high-k `k15` | Informative but still mixed | `Pantala flavescens` in both masked and unmasked matrices; `Campodea` closer than `Catajapyx` | Low-to-moderate; masked and unmasked broadly agree on the local `protura` neighborhood | No | No clean support; local neighborhood includes `Campodea`, but not a direct Diplura-only signal | Better than small-k, but still not taxonomically clean | Yes, as the current best diagnostic signal class, with caution |
| Phase 2 high-k `k21` | Saturated | Nearest-neighbor calls vary and are not reliable | Differences exist, but under saturation | No interpretable support | No interpretable support | No interpretable coherence | No; retain only as a saturation control |

## Main Conclusion

Conservative reading from existing outputs:

- Phase 2 small-k does not support a robust `Protura + Collembola` signal.
- Phase 2 high-k `k15` does not preserve the Phase 1 `protura_purged + Catajapyx` signal.
- Phase 2 high-k `k21` remains saturated and should remain declassed.
- High-k `k15` is more informative than small-k because it escapes the universal `protura_purged -> Daphnia` attraction seen in Phase 2 small-k.
- Even so, high-k `k15` is still not taxonomically clean: it places `protura_purged` in a mixed local neighborhood involving `Pantala`, `Campodea`, and `Thermobia`, not a clean diagnostic clade.
- At the current state, the workflow is strongest as a diagnostic framework for signal behavior and sensitivity analysis, not yet as a tree-first method that can support strong biological interpretation.

## Signal-by-Signal Comparison

### Phase 2 small-k

Using [phase2_smallk_diagnostics.md](/media/shared1/structural-phylogenomics-workflow/docs/datasets/phase2_smallk_diagnostics.md) plus the Phase 2 tree reports:

- `protura_purged` is nearest to `Daphnia pulex` in all `16/16` distance matrices.
- There is no stable `protura_purged + Folsomia candida` tree pattern.
- There is no stable `protura_purged + Diplura` pattern.
- Only one tree, `unmasked k8 cosine`, gives a direct `protura_purged + Campodea augens` cherry.
- `masked` trees are more stable than `unmasked` trees.
- `unmasked cosine` is the least stable small-k branch.

Interpretation:

- expanded taxon sampling breaks the apparent Phase 1 small-k `Protura + Folsomia` signal
- the dominant small-k effect is a likely non-phylogenetic attraction toward `Daphnia`

### Phase 2 high-k `k15`

Using [phase2_highk_diagnostics.md](/media/shared1/structural-phylogenomics-workflow/docs/datasets/phase2_highk_diagnostics.md) plus existing sketch outputs:

- `protura_purged` is nearest to `Pantala flavescens` in both masked and unmasked `k15` matrices.
- `Campodea augens` is closer to `protura_purged` than `Catajapyx aquilonaris` in both masked and unmasked `k15` matrices.
- `Catajapyx` is not the nearest neighbor of `protura_purged`.
- `protura_purged` does not form a direct cherry with either `Campodea` or `Catajapyx`.
- `protura_purged` does not form a direct cherry with any collembolan.
- Masked and unmasked `k15` agree on the broad local structure: `protura_purged` falls into a mixed local neighborhood containing `Campodea`, `Pantala`, and `Thermobia`.

Interpretation:

- Phase 1 `protura_purged + Catajapyx` does not survive expanded taxon sampling
- the current high-k signal is more informative than small-k, but still mixed and not taxonomically decisive

### Phase 2 high-k `k21`

Using the same existing high-k outputs:

- all pairwise distances are `>= 0.99` in both masked and unmasked `k21`
- mean distances are about `0.998-0.999`
- `k21` nearest-neighbor calls vary by dataset
- `k21` cherries differ between masked and unmasked trees

Interpretation:

- `k21` remains a saturation control only
- it should not be used to choose biological hypotheses or guide method development beyond documenting saturation

## Outlier Hypotheses Based Only on Existing Outputs

These are candidate drivers suggested by existing files only. They are not yet tested by dedicated sensitivity runs.

### Pantala flavescens

- In Phase 2 high-k `k15`, `Pantala flavescens` is the nearest neighbor of `protura_purged` in both masked and unmasked matrices.
- This makes `Pantala` the clearest candidate for leave-one-out testing.

### Thermobia domestica

- `Thermobia domestica` is extremely large.
- In high-k `k15`, it sits in the immediate local neighborhood around `protura_purged`.
- This makes it a plausible driver of mixed high-k local structure.

### Campodea augens

- `Campodea` is closer to `protura_purged` than `Catajapyx` in both high-k `k15` matrices.
- But this does not become a clean Diplura-wide signal.
- It may be acting as a lineage-specific or assembly-specific attractor rather than as evidence for a stable `Protura + Diplura` pattern.

### Catajapyx aquilonaris

- `Catajapyx` was the Phase 1 high-k nearest-neighbor anchor.
- In Phase 2 high-k `k15`, it loses that role.
- This makes it a direct taxon-sampling sensitivity case.

### Daphnia pulex

- `Daphnia` dominates Phase 2 small-k nearest-neighbor behavior in `16/16` matrices.
- That makes it the strongest candidate driver of small-k instability or compositional attraction.

### Podura aquatica

- `Podura` has a high masking fraction.
- It is a reasonable candidate for sensitivity testing because masking-related behavior could interact with both small-k and high-k signal structure.

### Fragmentation and genome-size heterogeneity

- `protura_purged`, `Catajapyx`, and `Campodea` are fragmented
- `Thermobia` is huge
- panel-wide genome-size heterogeneity is substantial

These factors remain plausible technical drivers of the mixed signal currently observed.

## What Can Be Concluded from Existing Files

### Supported by existing outputs

- Phase 2 small-k does not replicate a robust `protura_purged + Folsomia` signal.
- Phase 2 small-k is dominated by a stable `protura_purged -> Daphnia pulex` nearest-neighbor pattern.
- Phase 2 high-k `k15` does not replicate the Phase 1 `protura_purged + Catajapyx` signal.
- In Phase 2 high-k `k15`, `Campodea` is closer to `protura_purged` than `Catajapyx`, but neither yields a clean direct `Protura + Diplura` signal.
- Phase 2 high-k `k15` is more informative than Phase 2 small-k.
- Phase 2 high-k `k21` remains saturated and should remain declassed.
- The present workflow is better at diagnosing signal instability and taxon sensitivity than at yielding a trustworthy biological tree from this panel alone.

### Requires new jackknife or sensitivity runs

- Whether `Pantala flavescens` is the principal driver of the Phase 2 high-k `k15` local neighborhood
- Whether `Thermobia domestica` is distorting high-k neighborhood structure because of genome size
- Whether `Daphnia pulex` is the principal driver of the Phase 2 small-k pattern
- Whether `Podura aquatica` masking behavior contributes materially to the observed topology differences
- Whether removing `Catajapyx` or `Campodea` changes the apparent high-k Protura-adjacent signal
- Whether any single taxon or small taxon subset is responsible for most of the Phase 2 discordance

## Recommended Next Experiment

The next experiment should be a targeted taxon jackknife or outlier-sensitivity step using the current Phase 2 panel.

Recommended future leave-one-out candidates:

- without `Pantala flavescens`
- without `Thermobia domestica`
- without `Daphnia pulex`
- without `Podura aquatica`
- without `Catajapyx aquilonaris`
- without `Campodea augens`

Purpose of that step:

- test whether the current high-k `k15` neighborhood is driven by one or two influential taxa
- test whether the current small-k `Protura -> Daphnia` pattern collapses once likely outliers are removed
- determine whether a more taxonomically coherent signal emerges before adding new methodological layers

This should be treated as a technical sensitivity analysis, not as a biological adjudication.

## What Not To Do Yet

Based on current files, these should remain deferred:

- FASTQ
- simulated reads
- RepeatMasker or EDTA integration
- frequency filtering
- `k17`
- sketch tuning
- strong biological claims

The current bottleneck is not missing feature breadth. It is unresolved signal stability under the existing expanded assembly panel.

## Bottom Line

Phase 2 resolves one question clearly:

- neither the Phase 1 small-k `Protura + Folsomia` pattern nor the Phase 1 high-k `Protura + Catajapyx` pattern survives cleanly under expanded taxon sampling

What remains is a useful diagnostic outcome:

- small-k is unstable and likely compositionally distorted
- high-k `k15` is the best current signal class, but still taxonomically mixed
- `k21` is saturated and should stay declassed

The correct next step is not broader feature expansion. It is a controlled Phase 2 taxon jackknife/outlier sensitivity analysis using the workflow as it already exists.
