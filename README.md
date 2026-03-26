# Structural Phylogenomics Workflow

Snakemake workflow for an accession-driven, alignment-free phylogenomics pipeline on eukaryotic genomes. The current implementation covers metadata resolution, genome download, preprocessing, repeat masking, small-k spectrum generation, distance calculation, tree inference and tree comparison.

The default repeat backend uses `dustmasker` as a lightweight MVP. This is useful to get the workflow running, but it is not a substitute for a curated `RepeatModeler/RepeatMasker` TE annotation workflow.

## Goal

The project is designed to start from a minimal user input, a list of assembly accession IDs, and produce standardized metadata tables plus downstream phylogenomic outputs:

- cleaned and masked genomes
- k-mer feature matrices
- distance matrices
- Newick trees
- comparison tables across datasets and parameter combinations

## Input

Primary input:

- `metadata/accessions.txt`

One accession per line, for example:

```text
GCA_000934665.2
GCA_021134715.1
GCF_943734735.2
GCA_052721305.1
```

Alternative input from the command line:

```bash
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 4 --config accessions="GCA_000934665.2,GCA_021134715.1"
```

## Workflow Steps

### 1. Metadata Resolution

Input:

- accession list
- NCBI assembly summary tables

Actions:

- resolve each accession against NCBI assembly metadata
- collect assembly-level metadata
- group assemblies at organism level
- build a download manifest

Main outputs:

- `results/metadata/assemblies.tsv`
- `results/metadata/organisms.tsv`
- `results/metadata/download_manifest.tsv`

### 2. Genome Download

Input:

- `results/metadata/download_manifest.tsv`

Actions:

- download assembly FASTA files from NCBI FTP
- normalize local filenames as accession-based paths

Main outputs:

- `data/genomes/{accession}.fna.gz`

### 3. Raw Genome QC

Input:

- downloaded genomes

Actions:

- compute sequence count
- compute total assembly length
- compute length distribution and N50
- compute GC and `N` fraction

Main outputs:

- `results/qc/{accession}.tsv`
- `results/qc/qc_summary.tsv`

### 4. Genome Preprocessing

Input:

- downloaded genomes

Actions:

- remove short contigs
- optionally normalize headers
- optionally uppercase sequences
- optionally exclude organelle-like sequences by keyword

Main outputs:

- `results/preprocessed/{accession}.fa`

### 5. Repeat Annotation

Input:

- preprocessed genomes

Actions:

- annotate low-complexity / repeat-sensitive intervals with `dustmasker`

Main outputs:

- `results/repeats/annotation/{accession}.intervals.txt`
- `results/repeats/annotation/{accession}.summary.tsv`
- `results/repeats/repeat_annotation_summary.tsv`

### 6. Repeat Masking

Input:

- preprocessed genomes

Actions:

- generate masked genome FASTA files

Main outputs:

- `results/repeats/masked/{accession}.fa`

### 7. K-mer Spectrum Generation

Input:

- `unmasked` dataset: preprocessed genomes
- `masked` dataset: repeat-masked genomes

Actions:

- compute canonical k-mer spectra for configurable `k`
- export per-sample spectra
- assemble feature matrices per dataset and per `k`

Main outputs:

- `results/kmers/spectra/{dataset}/k{k}/{accession}.tsv`
- `results/kmers/matrices/{dataset}/k{k}.tsv`

### 8. Distance Matrix Generation

Input:

- k-mer feature matrices

Actions:

- compute pairwise distances between genomes
- currently supported metrics:
  - `cosine`
  - `jensen_shannon`

Main outputs:

- `results/distances/{dataset}/k{k}/{metric}.tsv`

### 9. Tree Inference

Input:

- distance matrices

Actions:

- infer trees from each matrix
- currently supported methods:
  - `nj`
  - `upgma`

Main outputs:

- `results/trees/{dataset}/k{k}/{metric}/{method}.nwk`

### 10. Tree Summary and Comparison

Input:

- all inferred trees

Actions:

- index produced trees
- compare tree pairs using internal split overlap
- report split-set Jaccard similarity

Main outputs:

- `results/reports/tree_manifest.tsv`
- `results/reports/tree_comparisons.tsv`

## Configuration

Main configuration file:

- `config/config.yaml`

Current configurable sections:

- `metadata`
- `downloads`
- `preprocessing`
- `repeat_annotation`
- `kmers`
- `distances`
- `trees`

## Repository Layout

- `metadata/`: user-provided accession lists
- `data/genomes/`: downloaded assemblies
- `config/`: workflow configuration
- `workflow/rules/`: Snakemake rules by module
- `workflow/scripts/`: Python helper scripts
- `resources/`: downloaded reference tables and local resources
- `results/`: generated workflow outputs

## Run

```bash
cd /home/bosi/kmer_phylo_workflow
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake -n
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 4
```

## Current Status

Implemented:

- accession-driven metadata resolution
- genome download
- QC
- preprocessing
- repeat annotation and masking
- small-k spectra
- distance matrices
- tree inference
- tree comparison reports

Not implemented yet:

- bootstrap / jackknife
- sliding-window phylogenomics
- high-k sketch-based comparisons
- tree visualization figures
- full TE-aware repeat annotation backend

## Notes

- Samples are keyed by assembly accession.
- Metadata are resolved from NCBI assembly summary tables.
- Downloaded genome filenames are normalized as `data/genomes/{accession}.fna.gz`.
- The current machine has a broken `/usr/local/bin/snakemake`; use `/home/bosi/miniforge3/envs/ampwrap/bin/snakemake`.

## Repository Scope

Track code, configuration and lightweight examples in Git.

Do not track:

- workflow outputs in `results/`
- Snakemake runtime state in `.snakemake/`
- downloaded genomes in `data/genomes/`
- large local resources in `resources/`
