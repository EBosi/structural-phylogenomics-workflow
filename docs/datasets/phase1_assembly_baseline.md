# Phase 1 Assembly Baseline Dataset

## Purpose

This dataset is the controlled assembly-based baseline for Phase 1 of the
structural phylogenomics roadmap. Its purpose is to test whether the existing
assembly workflow can produce interpretable QC, preprocessing, organelle
filtering, low-complexity masking, k-mer matrices, distance matrices and
assembly-based trees on a small real dataset before adding read-based analysis
or FASTQ support.

The dataset is intended for technical and preliminary biological diagnostics
only. It is not a formal reference phylogeny.

## Tracked Configuration

- Baseline config override: `config/phase1_assembly_baseline.yaml`
- Baseline dataset notes: `docs/datasets/phase1_assembly_baseline.md`

The config override only changes the local genome metadata table used for this
baseline. It does not change k-mer values, distance metrics, tree methods,
repeat masking settings or other biological parameters from `config/config.yaml`.

## Local Non-Tracked Paths

These paths are intentionally not tracked in Git:

- Source/staged assembly archive: `resources/phase1_assembly_baseline/assemblies/`
- Baseline provenance table: `resources/phase1_assembly_baseline/assembly_sources.tsv`
- Baseline local genome metadata table: `resources/phase1_assembly_baseline/local_genomes.tsv`
- Workflow-staged genomes: `data/genomes/`
- Workflow outputs: `results/`
- Organelle reference cache and BLAST database: `resources/organelle/`

The original local target genome used to create the staged assembly archive was:

- `/home/ebosi/shared1/protura_purged.fna.gz`

## Samples

| Sample/accession ID | Organism name | Biological role | Source/provenance | Known issues |
| --- | --- | --- | --- | --- |
| `protura_purged` | `protura_purged` | target | local assembly from `/home/ebosi/shared1/protura_purged.fna.gz` | Organism name/taxonomy to be confirmed; large genome; 8129 contigs; N50 186,972; no confident organelle contigs removed using the available mitochondrial references; highest low-complexity masked fraction in this dataset. |
| `GCA_000934665.2` | `Catajapyx aquilonaris (hexapods)` | close comparator, Diplura-like/entognath role to be confirmed | NCBI GenBank, `GCA_000934665.2_Caqu_2.0_genomic.fna.gz` | Scaffold-level assembly; highly fragmented; 20,583 raw sequences; 15,054 retained after min-length preprocessing; N50 45,890; no confident organelle contigs removed. |
| `GCA_052721305.1` | `Folsomia candida (springtails)` | close comparator, Collembola | NCBI GenBank, `GCA_052721305.1_ASM5272130v1_genomic.fna.gz` | Complete Genome assembly; no confident organelle contigs removed; role relative to target to be confirmed with reference tree. |
| `GCA_021134715.1` | `Daphnia pulex (common water flea)` | distant outgroup/control | NCBI RefSeq FASTA resolved from requested accession to `GCF_021134715.1_ASM2113471v1_genomic.fna.gz` | Requested accession was `GCA_021134715.1`; NCBI source FASTA used the RefSeq path `GCF_021134715.1`; one 15,333 bp confident mitochondrial contig removed. |
| `GCF_943734735.2` | `Anopheles gambiae (African malaria mosquito)` | outgroup/control, insect/hexapod | NCBI RefSeq, `GCF_943734735.2_idAnoGambNW_F1_1_genomic.fna.gz` | One 15,363 bp confident mitochondrial contig removed; biological role as outgroup/control should be formalized with a reference tree. |

## Commands Used For `pre_kmer`

The initial all-local-directory dry-run was:

```bash
snakemake -n --config local_genomes_dir=resources/phase1_assembly_baseline/assemblies
```

That command produced a valid DAG, but the corresponding `pre_kmer` run failed
at mitochondrial BLAST database construction because `local_genomes_dir` assigns
sample IDs as organism names. The organelle reference fetch therefore queried
organisms such as `GCA_000934665.2` instead of real organism names and produced
an empty mitochondrial reference FASTA.

The successful `pre_kmer` run used a local metadata TSV with organism names:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 pre_kmer \
  --config metadata='{"accession_file":"metadata/accessions.txt","local_genomes_file":"resources/phase1_assembly_baseline/local_genomes.tsv","local_genomes_dir":"","local_genome_extensions":[".fa",".fasta",".fna",".fa.gz",".fasta.gz",".fna.gz"],"ncbi_eutils_base":"https://eutils.ncbi.nlm.nih.gov/entrez/eutils","request_timeout_seconds":60}'
```

For a fresh run, use the tracked config override instead of the inline JSON:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 pre_kmer \
  --configfile config/config.yaml \
  --configfile config/phase1_assembly_baseline.yaml
```

When continuing from the already completed baseline outputs, use mtime-based
rerun triggers to avoid recomputing `pre_kmer` only because the original run used
inline JSON and Snakemake recorded a different provenance representation for
`metadata.request_timeout_seconds`:

```bash
conda run -n structural-phylogenomics-workflow snakemake --cores 4 \
  results/reports/tree_manifest.tsv \
  results/reports/tree_comparisons.tsv \
  --configfile config/config.yaml \
  --configfile config/phase1_assembly_baseline.yaml \
  --rerun-triggers mtime
```

## Main Outputs Generated

Pre-kmer outputs:

- `results/metadata/assemblies.tsv`
- `results/metadata/organisms.tsv`
- `results/metadata/download_manifest.tsv`
- `results/qc/qc_summary.tsv`
- `results/preprocessing/preprocessing_summary.tsv`
- `results/organelle/organelle_summary.tsv`
- `results/repeats/repeat_annotation_summary.tsv`
- `results/reports/pre_kmer_summary.tsv`
- `results/reports/pre_kmer_report.md`

Assembly-tree continuation outputs:

- `results/kmers/spectra/`
- `results/kmers/matrices/`
- `results/distances/`
- `results/trees/`
- `results/reports/tree_manifest.tsv`
- `results/reports/tree_comparisons.tsv`

The assembly-tree continuation used the default configured values from
`config/config.yaml`:

- datasets: `unmasked`, `masked`
- k values: `5`, `6`, `7`, `8`
- distance metrics: `cosine`, `jensen_shannon`
- tree methods: `nj`, `upgma`

## Baseline `pre_kmer` Summary

- Samples: 5
- Total raw bases: 1,709,478,007
- Total bases after preprocessing: 1,706,768,393
- Total bases after organelle filtering: 1,706,737,697
- Total bases removed as organellar: 30,696
- Low-complexity/repeat-sensitive masking backend: `dustmasker`

The dataset is usable with caution.

## What This Dataset Can Support

- Technical validation of the assembly-based workflow on a small real dataset.
- Preliminary checks for sample-level QC outliers.
- Preliminary comparison of assembly-based k-mer trees across configured
  datasets, k values, metrics and tree methods.
- Identification of obvious instability before adding a formal reference tree.

## What This Dataset Cannot Support

- Strong phylogenetic claims.
- Interpretation of discordant trees as biological signal.
- Replacement of a marker-based, BUSCO/single-copy ortholog or literature
  reference tree.
- Claims about read-based k-mer performance.
- Claims about FASTQ support.
- Claims that low-complexity masking is equivalent to curated TE annotation.
- Confident conclusions about organelle effects, because mitochondrial
  references were only recovered for some organisms in the baseline.
