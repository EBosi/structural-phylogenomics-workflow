# Phase 2 Expanded Assembly Panel

## Purpose

Phase 2 tests whether the Phase 1 discordance persists after adding broader assembly-based taxon sampling. Phase 1 used a small 5-taxon panel where small-k NJ mostly supported a `protura_purged + Folsomia` diagnostic hypothesis, high-k k15 supported a `protura_purged + Catajapyx` diagnostic hypothesis, and high-k k21 was declassed because distances were saturated.

This Phase 2 panel is diagnostic only. It is not a definitive phylogenetic dataset, and k-mer/sketch concordance with any topology must not be interpreted as biological proof without an external curated reference framework.

## Input and Output Directories

- Phase 2 resources: `resources/phase2_expanded_panel/`
- Phase 2 staged assemblies: `resources/phase2_expanded_panel/assemblies/`
- Phase 2 output root: `results/phase2_expanded_panel/`
- Phase 2 local genome table: `resources/phase2_expanded_panel/local_genomes.tsv`
- Phase 2 assembly source table: `resources/phase2_expanded_panel/assembly_sources.tsv`
- Phase 2 config override: `config/phase2_expanded_panel.yaml`

Transcriptomes are not mixed with genome assemblies in the main run. Protura and Diplura transcriptomes remain side-panel entries only.

## Why Phase 1 Was Insufficient

The Phase 1 panel had only one local Protura target, one Collembola comparator, one Diplura comparator, one insect comparator, and one non-hexapod arthropod outgroup. That sampling was too small to separate possible biological signal from compositional bias, assembly-specific artifacts, repeat content, genome size effects, or instability caused by small-k distances.

Phase 2 expands Collembola, Diplura, ectognath/insect comparators, and crustacean outgroups while keeping the analysis assembly-only.

## Selection Rules

Selected main-run taxa were taken from `docs/datasets/phase2_revised_panel.tsv` using:

- `include_phase2 == yes`
- `datatype` equal to `genome` or `genome_local`
- `decision == main_panel`, plus `decision == to_add` with `panel_tier == main_genome_tentative`

Excluded from the main run:

- transcriptomes
- reserve taxa
- `exclude_now` taxa

## Selected Main Panel

Main panel size: 16 genome assemblies.

Group counts:

- protura: 1
- collembola: 6
- diplura: 2
- insecta: 4
- crustacea: 3

| Sample ID | Source accession | Organism | Broad group | Priority | Role | Caveats |
| --- | --- | --- | --- | --- | --- | --- |
| protura_purged | local:protura_purged | protura_purged | protura | high | target | Local assembly; taxonomy and assembly quality must remain explicitly caveated |
| GCA_947179485.1 | GCA_947179485.1 | Allacma fusca | collembola | high | collembola_comparator | Check assembly size/fragmentation before final run |
| GCA_964034955.1 | GCA_964034955.1 | Orchesella flavescens | collembola | high | collembola_comparator | Check assembly size/fragmentation before final run |
| GCA_052721305.1 | GCA_052721305.1 | Folsomia candida | collembola | high | phase1_anchor_collembola | Keep for continuity with Phase 1 |
| GCA_046254475.1 | GCA_046254475.1 | Podura aquatica | collembola | medium | collembola_comparator | Medium priority due to BUSCO and possible fragmentation |
| GCA_965653795.1 | GCA_965653795.1 | Pogonognathellus longicornis | collembola | high | collembola_comparator | Check assembly size/fragmentation before final run |
| GCA_964276635.1 | GCA_964276635.1 | Lepidocyrtus curvicollis | collembola | high | collembola_comparator | Check assembly size/fragmentation before final run |
| GCA_000934665.2 | GCA_000934665.2 | Catajapyx aquilonaris | diplura | high | phase1_anchor_diplura | Single Diplura genome is not enough; pair with Campodea if possible |
| GCA_009757345.1 | GCA_009757345.1 | Campodea augens | diplura | high | second_diplura_candidate | Verified to NCBI GenBank FTP and included; BUSCO/quality still to verify |
| GCA_964235325.1 | GCA_964235325.1 | Thermobia domestica | insecta | high | ectognath_comparator | Check assembly size/fragmentation before final run |
| GCA_949628265.1 | GCF_949628265.1 | Cloeon dipterum | insecta | high | ectognath_comparator | Source FTP uses GCF_949628265.1 while the sample ID preserves the table accession style |
| GCA_020796165.1 | GCA_020796165.1 | Pantala flavescens | insecta | high | ectognath_comparator | Check assembly size/fragmentation before final run |
| GCF_943734735.2 | GCF_943734735.2 | Anopheles gambiae | insecta | high | phase1_anchor_insect | Use as ectognath/insect comparator, not strict Hexapoda outgroup |
| GCA_007210705.1 | GCF_007210705.1 | Tigriopus californicus | crustacea | high | outgroup_crustacean | Source FTP uses GCF_007210705.1 while the sample ID preserves the table accession style |
| GCA_021134715.1 | GCF_021134715.1 | Daphnia pulex | crustacea | high | phase1_anchor_outgroup | Source FTP uses GCF_021134715.1 while the sample ID preserves the table accession style |
| GCA_965194735.1 | GCA_965194735.1 | Jaera ischiosetosa | crustacea | medium | outgroup_crustacean | Medium priority; check assembly size/fragmentation before final run |

## Side, Reserve, and Excluded Taxa

Side transcriptome count: 6.

- GKDB00000000.1, Sinentomon erythranum, protura
- GAXE00000000.2, Acerentomon sp., protura
- GAXJ00000000.2, Occasjapyx japonicus, diplura
- GKDA00000000.1, Lepidocampa weberi, diplura
- GKCZ00000000.1, Octostigma sinensis, diplura
- HBDP00000000.1, Machilis pallida, insecta

Reserve genome count: 8.

- GCA_009870155.1, Pseudobourletiella spinata, collembola
- GCA_046579375.1, Orthonychiurus folsomi, collembola
- GCA_046563655.1, Yuukianura sp., collembola
- GCA_009869795.1, Neelides sp., collembola
- GCA_977018435.1, Lepisma saccharinum, insecta
- GCA_977061575.1, Ephemera vulgata, insecta
- GCA_977012765.1, Anax junius, insecta
- GCA_032884065.1, Artemia franciscana, crustacea

Excluded or deprioritized count: 9.

- GCA_054131365.1, Gryllus bimaculatus, insecta
- GCA_976940595.1, Empusa pennicornis, insecta
- GCA_050494535.1, Timema cristinae, insecta
- GCA_918797505.1, Bemisia tabaci, insecta
- GCA_053478235.1, Coccinella transversoguttata, insecta
- GCA_912999745.1, Papilio machaon, insecta
- GCA_910591885.2, Bombus terrestris, insecta
- GCA_053477335.1, Acartia tonsa, crustacea
- GCA_052075045.1, Palinurus elephas, crustacea

## Staging Status

All 16 main-panel assemblies are staged under `resources/phase2_expanded_panel/assemblies/`.

- `protura_purged` was reused by symlink, not copied.
- Phase 1 anchors `GCA_052721305.1`, `GCA_000934665.2`, `GCF_943734735.2`, and `GCA_021134715.1` were reused by symlink.
- Eleven public assemblies were downloaded directly under the Phase 2 assembly directory.
- `Campodea augens` / `GCA_009757345.1` was verified to a GenBank FASTA URL and staged.
- Gzip FASTA validation succeeded for the 16 staged assembly paths.

Missing or `to_verify` taxa for the operational main run: none. Some assembly-level caveats remain to verify, especially BUSCO/completeness, fragmentation, and GCA/GCF accession consistency for entries where the source FTP uses GCF while the stable sample ID preserves the panel table ID.

## Config

The Phase 2 override is `config/phase2_expanded_panel.yaml`.

It sets:

- `output_root: "results/phase2_expanded_panel"`
- `genome_root: "resources/phase2_expanded_panel/assemblies"`
- `metadata.local_genomes_file: "resources/phase2_expanded_panel/local_genomes.tsv"`
- `metadata.local_genomes_dir: ""`
- `trees.methods: ["nj"]`
- `sketch.datasets: ["masked", "unmasked"]`
- `sketch.k_values: [15, 17, 21]`
- `sketch.methods: ["nj"]`

The main `config/config.yaml` was not modified.

## Commands Run

Dry-run for pre-kmer summary:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/phase2_expanded_panel/reports/pre_kmer_summary.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

Dry-run for high-k sketch k15 masked NJ:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/phase2_expanded_panel/sketch/trees/masked/k15/nj.nwk \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

Dry-run for small-k tree reports:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/phase2_expanded_panel/reports/tree_manifest.tsv \
  results/phase2_expanded_panel/reports/tree_comparisons.tsv \
  --configfile config/config.yaml \
  --configfile config/phase2_expanded_panel.yaml \
  --rerun-triggers mtime
```

## Dry-Run Result and Current Blocker

The dry-runs showed that the requested analysis outputs under `results/` would be rooted under `results/phase2_expanded_panel/`, and genome inputs came from `resources/phase2_expanded_panel/assemblies/`.

However, the DAG still plans global organelle resource outputs:

- `resources/organelle/mitochondrion_refs.fasta`
- `resources/organelle/mitochondrion_refs.tsv`
- `resources/organelle/blastdb/mitochondrion.nhr`
- `resources/organelle/blastdb/mitochondrion.nin`
- `resources/organelle/blastdb/mitochondrion.nsq`

Because the Phase 2 guardrail requires dedicated Phase 2 input/output directories, no real Phase 2 analysis was run. The next infrastructure decision is whether to add a minimal configurable organelle/resource root, reuse pre-existing immutable organelle references explicitly, or otherwise isolate these resource outputs before running Phase 2.

## Analysis Status

No real Phase 2 k-mer, distance, tree, or sketch analysis has been executed yet.

Therefore:

- pre-kmer quality summary: not produced
- small-k result summary: not produced
- high-k k15/k17/k21 saturation summary: not produced
- small-k versus high-k Phase 2 agreement/discordance: not assessed

This is intentional. Running the analysis before isolating the global organelle resource outputs would risk mixing Phase 2 state with existing global resources.

## Caveats

- The expanded panel reduces, but does not eliminate, taxon-sampling limitations.
- Public assembly quality, completeness, fragmentation, haplotig/purge status, contamination, and repeat content remain potential confounders.
- `protura_purged` taxonomy and assembly quality remain to be confirmed.
- `Anopheles gambiae` is an ectognath/insect comparator, not a strict Hexapoda outgroup.
- `Daphnia pulex`, `Tigriopus californicus`, and `Jaera ischiosetosa` are non-hexapod arthropod/crustacean outgroups.
- No strong biological claims should be made from Phase 2 until the isolated run completes and results are compared against a curated biological reference framework.
