# Phase 1 High-k Sketch Diagnostics

## Scope

This document records the approved Phase 1 high-k sketch micro-experiment.
The goal was to test whether masked assembly MinHash/Jaccard at higher k
reduces the unstable compositional signal observed in the small-k k-mer trees.

This run did not modify workflow code, did not change `config/config.yaml`, did
not run `full_analysis`, and did not introduce FASTQ or simulated reads.

Interpretation is conservative: these outputs are a diagnostic comparison
against the small-k baseline, not conclusive phylogenetic evidence.

## Config Used

Dedicated override:

`config/phase1_highk_sketch.yaml`

```yaml
metadata:
  local_genomes_file: "resources/phase1_assembly_baseline/local_genomes.tsv"
  local_genomes_dir: ""

sketch:
  datasets: ["masked"]
  k_values: [15, 21]
  methods: ["nj"]
```

Default values retained from `config/config.yaml`:

| Parameter | Value |
| --- | --- |
| sketch input dataset | `masked` |
| k values | `15`, `21` |
| tree method | `nj` |
| distance | `minhash_jaccard` |
| canonical k-mers | `true` |
| requested MinHash hashes | `1000` |

## Commands Executed

Dry-run:

```bash
conda run -n structural-phylogenomics-workflow snakemake -n \
  results/sketch/trees/masked/k15/nj.nwk \
  results/sketch/trees/masked/k21/nj.nwk \
  results/sketch/distances/masked/k15/minhash_jaccard.tsv \
  results/sketch/distances/masked/k21/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile config/phase1_assembly_baseline.yaml \
  --configfile config/phase1_highk_sketch.yaml \
  --rerun-triggers mtime
```

Real run:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 \
  results/sketch/trees/masked/k15/nj.nwk \
  results/sketch/trees/masked/k21/nj.nwk \
  results/sketch/distances/masked/k15/minhash_jaccard.tsv \
  results/sketch/distances/masked/k21/minhash_jaccard.tsv \
  --configfile config/config.yaml \
  --configfile config/phase1_assembly_baseline.yaml \
  --configfile config/phase1_highk_sketch.yaml \
  --rerun-triggers mtime
```

The real run completed 14 Snakemake jobs:

| Rule | Count |
| --- | ---: |
| `sketch_signature_per_sample` | 10 |
| `sketch_distance_matrix` | 2 |
| `infer_sketch_tree` | 2 |

## Outputs Produced

Signature files:

| File | Effective hashes | Size bytes |
| --- | ---: | ---: |
| `results/sketch/signatures/masked/k15/GCA_000934665.2.tsv` | 1000 | 46062 |
| `results/sketch/signatures/masked/k15/GCA_021134715.1.tsv` | 1000 | 46448 |
| `results/sketch/signatures/masked/k15/GCA_052721305.1.tsv` | 1000 | 46226 |
| `results/sketch/signatures/masked/k15/GCF_943734735.2.tsv` | 1000 | 46081 |
| `results/sketch/signatures/masked/k15/protura_purged.tsv` | 1000 | 44958 |
| `results/sketch/signatures/masked/k21/GCA_000934665.2.tsv` | 1000 | 45804 |
| `results/sketch/signatures/masked/k21/GCA_021134715.1.tsv` | 1000 | 46283 |
| `results/sketch/signatures/masked/k21/GCA_052721305.1.tsv` | 1000 | 45813 |
| `results/sketch/signatures/masked/k21/GCF_943734735.2.tsv` | 1000 | 45826 |
| `results/sketch/signatures/masked/k21/protura_purged.tsv` | 1000 | 44752 |

Distance matrices:

| File | Distance range | Mean pairwise distance |
| --- | --- | ---: |
| `results/sketch/distances/masked/k15/minhash_jaccard.tsv` | 0.742138-0.889506 | 0.817997 |
| `results/sketch/distances/masked/k21/minhash_jaccard.tsv` | 0.998999-1.000000 | 0.999800 |

Trees:

| File |
| --- |
| `results/sketch/trees/masked/k15/nj.nwk` |
| `results/sketch/trees/masked/k21/nj.nwk` |

## Nearest Neighbors

Nearest neighbors from `results/sketch/distances/masked/k15/minhash_jaccard.tsv`:

| Sample | Nearest neighbor | Distance |
| --- | --- | ---: |
| `protura_purged` | `GCA_000934665.2` | 0.742138 |
| `GCA_000934665.2` | `protura_purged` | 0.742138 |
| `GCA_021134715.1` | `protura_purged` | 0.853868 |
| `GCF_943734735.2` | `protura_purged` | 0.774510 |
| `GCA_052721305.1` | `protura_purged` | 0.783455 |

Nearest neighbors from `results/sketch/distances/masked/k21/minhash_jaccard.tsv`:

| Sample | Nearest neighbor | Distance |
| --- | --- | ---: |
| `protura_purged` | `GCA_000934665.2` | 1.000000 |
| `GCA_000934665.2` | `GCA_052721305.1` | 0.998999 |
| `GCA_021134715.1` | `protura_purged` | 1.000000 |
| `GCF_943734735.2` | `GCA_052721305.1` | 0.998999 |
| `GCA_052721305.1` | `GCA_000934665.2` | 0.998999 |

For k21, most distances are exactly `1.000000`; nearest-neighbor calls at this
k are therefore dominated by ties or near-ties and should not be interpreted as
stable biological signal.

## Topologies

k15 NJ tree:

```text
(GCA_052721305.1:0.39244,(GCF_943734735.2:0.39981,(protura_purged:0.35818,GCA_000934665.2:0.38396)Inner1:0.02043)Inner2:0.00757,GCA_021134715.1:0.46861)Inner3:0.00000;
```

Summary: k15 places `protura_purged` with `GCA_000934665.2`
(`Catajapyx aquilonaris`) as a cherry, with `GCF_943734735.2`
(`Anopheles gambiae`) attaching next in the NJ topology. `GCA_052721305.1`
(`Folsomia candida`) does not form a cherry with `protura_purged`.

k21 NJ tree:

```text
((GCA_021134715.1:0.50000,protura_purged:0.50000)Inner2:0.00025,(GCA_052721305.1:0.49933,GCA_000934665.2:0.49967)Inner1:0.00025,GCF_943734735.2:0.49975)Inner3:0.00000;
```

Summary: k21 places `GCA_021134715.1` (`Daphnia pulex`) with
`protura_purged`, and `GCA_052721305.1` (`Folsomia candida`) with
`GCA_000934665.2` (`Catajapyx aquilonaris`). Because the distance matrix is
nearly saturated, this topology is not reliable.

## Comparison With Small-k Baseline

Question: does `Daphnia pulex` remain nearest neighbor of `protura_purged`?

Answer: no for k15; `protura_purged` is nearest to `GCA_000934665.2`
(`Catajapyx aquilonaris`). For k21 all `protura_purged` off-diagonal distances
are `1.000000`, so the reported nearest neighbor is an order/tie artifact.

Question: does `Folsomia candida` + `protura_purged` emerge?

Answer: no. k15 does not place `Folsomia candida` with `protura_purged`, and k21
is saturated. This weakens the specific small-k NJ signal where
`Folsomia candida` + `protura_purged` appeared repeatedly.

Question: does `Catajapyx aquilonaris` continue to attach to
`Anopheles gambiae`?

Answer: not as the main cherry. In k15, `Catajapyx aquilonaris` is closest to
`protura_purged`, with `Anopheles gambiae` attaching outside that pair. In k21,
`Catajapyx aquilonaris` pairs with `Folsomia candida`, but the matrix is nearly
saturated and should not be overinterpreted.

Question: are k15 and k21 coherent?

Answer: no. k15 has a non-saturated distance range and gives one topology;
k21 is almost fully saturated and gives a different topology. The discordance is
most parsimoniously technical, caused by loss of shared high-k MinHash signal in
this small and phylogenetically broad assembly panel.

## Saturation Check

k15 distances are high but still variable:

| Statistic | Value |
| --- | ---: |
| minimum | 0.742138 |
| maximum | 0.889506 |
| mean | 0.817997 |

k21 distances are almost completely saturated:

| Statistic | Value |
| --- | ---: |
| minimum | 0.998999 |
| maximum | 1.000000 |
| mean | 0.999800 |

Interpretation: k15 retains some comparative information. k21 is not
informative with the current 1000-hash sketch and this five-assembly dataset.

## Judgment

Overall judgment: result ambiguous.

High-k sketch is technically functional and k15 reduces the specific small-k
artifact where `Daphnia pulex` was often nearest to `protura_purged`. However,
k15 does not recover the expected `Folsomia candida` + `protura_purged` signal,
and k21 is effectively saturated.

The result should not be treated as a strong phylogenetic result. It is useful
as a diagnostic showing that k21 MinHash/Jaccard at the current sketch size is
too sparse for this panel, while k15 may remain worth comparing against a
provisional reference tree.

The next appropriate step is to build and compare against a provisional
reference tree. Do not interpret discordances biologically until that comparison
exists.
