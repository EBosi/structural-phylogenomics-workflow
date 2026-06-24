# Phase 1 Taxonomic Reference Notes

## Scope

This note provides a cautious taxonomic framework for interpreting the Phase 1
diagnostic discordance:

- small-k NJ mostly supports Hypothesis A: `protura_purged` + Folsomia
- high-k k15 supports Hypothesis B: `protura_purged` + Catajapyx
- high-k k21 shows Hypothesis C, but k21 is saturated and declassed

This document does not introduce new analyses, new code, new dependencies or a
definitive reference phylogeny. It records what should be considered before
interpreting the current k-mer/sketch results biologically.

## Sample Table

| Sample ID | Organism name | Broad taxonomic group | Role in panel | Caveats |
| --- | --- | --- | --- | --- |
| `protura_purged` | unknown / to be confirmed | Protura candidate / target taxon to confirm | target assembly | Exact taxonomy and provenance need confirmation. Assembly-specific effects, genome size, fragmentation, repeat content and masking fraction may influence k-mer/sketch signal. |
| `GCA_052721305.1` | `Folsomia candida` | Collembola | Collembola comparator | Used as a close hexapod comparator for an Ellipura-like signal. Exact phylogenetic interpretation needs external curated support; source to verify. |
| `GCA_000934665.2` | `Catajapyx aquilonaris` | Diplura | Diplura comparator | Used as a Diplura comparator for a Nonoculata-like signal. The assembly is fragmented in this baseline, so concordance may reflect technical effects as well as biological signal. |
| `GCF_943734735.2` | `Anopheles gambiae` | Insecta / ectognath hexapod | ectognath comparator / distant hexapod control | Not a strict outgroup to Hexapoda. It is a hexapod insect comparator and should not be treated as an external arthropod outgroup. |
| `GCA_021134715.1` | `Daphnia pulex` | Crustacea / non-hexapod arthropod | non-hexapod arthropod outgroup/control | Broad arthropod outgroup relative to the hexapod-focused samples. Its behavior in small-k nearest-neighbor diagnostics was suspicious and may reflect compositional signal. |

## Clarified Taxonomic Roles

`protura_purged` is the target assembly. Its taxonomy is still to be confirmed.
All biological interpretation depends on validating this target identity and
checking whether the assembly is representative and clean enough for comparison.

`GCA_052721305.1` / `Folsomia candida` is the Collembola comparator. It is used
to test whether Phase 1 trees show an Ellipura-like signal, meaning a
Protura + Collembola affinity. Exact literature support and the expected
relationship in the specific taxonomic framework should be verified externally;
source to verify.

`GCA_000934665.2` / `Catajapyx aquilonaris` is the Diplura comparator. It is
used to test whether Phase 1 trees show a Nonoculata-like signal, meaning a
Protura + Diplura affinity. Exact literature support and competing hypotheses
should be verified externally; source to verify.

`GCF_943734735.2` / `Anopheles gambiae` is an insect / ectognath hexapod
comparator. It is useful as a distant hexapod control, but it is not a strict
outgroup to Hexapoda because it is itself a hexapod.

`GCA_021134715.1` / `Daphnia pulex` is a non-hexapod arthropod / crustacean
outgroup. It provides broader arthropod context, but the five-taxon panel is too
small to robustly place deep arthropod or hexapod relationships.

## Diagnostic Hypotheses

### A: Protura + Collembola / Ellipura-like Signal

Hypothesis A groups `protura_purged` with `GCA_052721305.1`
(`Folsomia candida`).

What it means diagnostically:

This is the hypothesis most consistent with the current small-k NJ results.
Small-k NJ trees contain the focal split `protura_purged` + Folsomia in 15/16
comparisons.

What it does not prove:

It does not prove an Ellipura relationship. Small-k spectra can be strongly
affected by base composition, low-order sequence composition, genome size,
repeat content, masking behavior and assembly properties.

Evidence status:

Diagnostic only. External taxonomic and phylogenomic sources need to be checked
before treating this as biological support; source to verify.

### B: Protura + Diplura / Nonoculata-like Signal

Hypothesis B groups `protura_purged` with `GCA_000934665.2`
(`Catajapyx aquilonaris`).

What it means diagnostically:

This is the hypothesis most consistent with the current high-k k15 sketch
results. Both masked and unmasked high-k k15 NJ trees contain the focal split
`protura_purged` + Catajapyx.

What it does not prove:

It does not prove a Nonoculata relationship. High-k sketch overlap may reflect
fragment sharing, assembly structure, repeat/low-complexity handling,
taxon-specific genome properties or other technical factors.

Evidence status:

Diagnostic only. External taxonomic and phylogenomic sources need to be checked
before treating this as biological support; source to verify.

### C: Collembola + Diplura / Protura-sister or Alternative Deep-Hexapod Signal

Hypothesis C groups `GCA_052721305.1` (`Folsomia candida`) with
`GCA_000934665.2` (`Catajapyx aquilonaris`).

What it means diagnostically:

This hypothesis would separate `protura_purged` from the two comparator
apterygote-like hexapod lineages in this small panel. It can be considered an
alternative deep-hexapod diagnostic pattern or a Protura-sister-style scenario,
depending on how an external reference is later defined.

What it does not prove:

The current apparent k21 support for C is not interpretable because the high-k
k21 distance matrix is nearly saturated. Saturated distances make nearest
neighbors and tree topology unstable or dominated by ties.

Evidence status:

Declassed for Phase 1 because the only current support comes from saturated
k21 sketch output. Any biological relevance requires external reference
assessment and broader taxon sampling; source to verify.

## Implications for Phase 1 Results

Small-k support for Hypothesis A may reflect a true Ellipura-like signal, but it
may also reflect compositional signal. The previous diagnostics showed that
small-k k-mer trees are unstable across parameters and likely sensitive to
low-order genome composition.

High-k k15 support for Hypothesis B may reflect a true Nonoculata-like signal,
but it may also reflect fragment-sharing signal, assembly structure, repeat
content or other high-k overlap effects. The masking sensitivity analysis showed
that masking changes high-k k15 distances while leaving the nearest-neighbor and
topology unchanged.

High-k k21 support for Hypothesis C is not interpretable here. The k21 MinHash
Jaccard distances are almost completely saturated, so the k21 tree should remain
declassed in Phase 1 interpretation.

The discordance between small-k A and high-k k15 B should not be treated as a
biological result yet. It is a diagnostic observation that motivates a curated
external reference framework and broader taxon sampling.

The next computational step should wait until this biological framework is
recorded and reviewed. In particular, FASTQ, simulated reads, serious repeat
masking and new workflow features should remain out of scope until there is a
clearer reference basis for interpreting the current assembly-based discordance.

## External Sources To Verify

The following points require external curated verification before they are used
as biological assumptions:

| Topic | Why it matters | Status |
| --- | --- | --- |
| Confirmed taxonomy and provenance of `protura_purged` | The target identity controls interpretation of every hypothesis. | source to verify |
| Current accepted placement of Protura relative to Collembola and Diplura | Needed to decide whether A or B is closer to an external reference. | source to verify |
| Evidence and terminology for Ellipura-like hypotheses | Needed before treating Protura + Collembola as a reference expectation. | source to verify |
| Evidence and terminology for Nonoculata-like hypotheses | Needed before treating Protura + Diplura as a reference expectation. | source to verify |
| Appropriate outgroup structure for a minimal hexapod/arthropod panel | Anopheles is not a strict Hexapoda outgroup; Daphnia is broader arthropod context. | source to verify |
| Whether additional taxa are needed before interpreting deep-hexapod signal | The current five-taxon panel is too small for strong phylogenetic claims. | source to verify |

## Conservative Conclusion

The Phase 1 reference hypotheses are useful for organizing the discordance, but
they are not a definitive reference tree. Concordance between a k-mer/sketch tree
and one hypothesis is not proof of biological truth.

At this stage, the most defensible conclusion is methodological: small-k NJ and
high-k k15 emphasize different signals. Interpreting that discordance requires
an external curated phylogeny or broader taxon sampling before stronger
biological claims are made.
