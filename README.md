# Structural Phylogenomics Workflow

Snakemake workflow for an accession-driven, alignment-free phylogenomics pipeline on eukaryotic genomes. The current implementation covers metadata resolution, genome download, QC, preprocessing, organelle screening, low-complexity masking, small-k spectrum generation, distance calculation, tree inference and tree comparison.

The default repeat backend uses `dustmasker` as a lightweight MVP. This is useful to get the workflow running, but it is not a substitute for a curated `RepeatModeler/RepeatMasker` TE annotation workflow.

## Goal

The project is designed to start from a minimal user input, a list of assembly accession IDs, and produce standardized metadata tables plus downstream phylogenomic outputs:

- cleaned, organelle-filtered and masked genomes
- k-mer feature matrices
- distance matrices
- Newick trees
- comparison tables across datasets and parameter combinations
- bootstrap and jackknife support summaries
- high-k MinHash sketch distances and trees

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
- NCBI E-utilities

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

Main outputs:

- `results/preprocessed/{accession}.fa`
- `results/preprocessing/{accession}.summary.tsv`
- `results/preprocessing/preprocessing_summary.tsv`

### 5. Organelle Screening and Filtering

Input:

- preprocessed genomes
- mitochondrial reference sequences fetched from NCBI

Actions:

- build a local BLAST database of mitochondrial references
- screen contigs against organelle references
- classify contigs as `nuclear_like`, `organelle_ambiguous` or `organelle_confident`
- remove only `organelle_confident` contigs

Main outputs:

- `results/organelle/calls/{accession}.tsv`
- `results/organelle/filtered/{accession}.fa`
- `results/organelle/{accession}.summary.tsv`
- `results/organelle/organelle_summary.tsv`

### 6. Low-Complexity Annotation

Input:

- organelle-filtered genomes

Actions:

- annotate low-complexity / repeat-sensitive intervals with `dustmasker`

Main outputs:

- `results/repeats/annotation/{accession}.intervals.txt`
- `results/repeats/annotation/{accession}.summary.tsv`
- `results/repeats/repeat_annotation_summary.tsv`

### 7. Low-Complexity Masking

Input:

- organelle-filtered genomes

Actions:

- generate masked genome FASTA files

Main outputs:

- `results/repeats/masked/{accession}.fa`

### 8. K-mer Spectrum Generation

Input:

- `unmasked` dataset: organelle-filtered genomes
- `masked` dataset: repeat-masked genomes

Actions:

- compute canonical k-mer spectra for configurable `k`
- export per-sample spectra
- assemble feature matrices per dataset and per `k`

Main outputs:

- `results/kmers/spectra/{dataset}/k{k}/{accession}.tsv`
- `results/kmers/matrices/{dataset}/k{k}.tsv`

### 9. Distance Matrix Generation

Input:

- k-mer feature matrices

Actions:

- compute pairwise distances between genomes
- currently supported metrics:
  - `cosine`
  - `jensen_shannon`

Main outputs:

- `results/distances/{dataset}/k{k}/{metric}.tsv`

### 10. Tree Inference

Input:

- distance matrices

Actions:

- infer trees from each matrix
- currently supported methods:
  - `nj`
  - `upgma`

Main outputs:

- `results/trees/{dataset}/k{k}/{metric}/{method}.nwk`

### 11. Tree Summary and Comparison

Input:

- all inferred trees

Actions:

- index produced trees
- compare tree pairs using internal split overlap
- report split-set Jaccard similarity

Main outputs:

- `results/reports/tree_manifest.tsv`
- `results/reports/tree_comparisons.tsv`

## Pre-kmer Milestone

The current default workflow target is the pre-kmer stage. This stage includes:

- metadata resolution
- genome download
- raw QC
- genome preprocessing
- organelle screening and filtering
- low-complexity annotation
- low-complexity masking
- consolidated pre-kmer reporting

Main milestone outputs:

- `results/reports/pre_kmer_summary.tsv`
- `results/reports/pre_kmer_report.md`

### 12. Window-Based Bootstrap and Contig Jackknife

Input:

- inferred reference trees
- partitioned k-mer spectra by genomic windows or contigs

Actions:

- compute window-level spectra for bootstrap
- compute contig-level spectra for jackknife
- resample units within each genome
- infer replicate trees
- estimate split support against the full-data reference tree

Main outputs:

- `results/resampling/bootstrap/{dataset}/k{k}/{metric}/{method}.support.tsv`
- `results/resampling/bootstrap/{dataset}/k{k}/{metric}/{method}.summary.tsv`
- `results/resampling/jackknife/{dataset}/k{k}/{metric}/{method}.support.tsv`
- `results/resampling/jackknife/{dataset}/k{k}/{metric}/{method}.summary.tsv`

Operational notes:

- run resampling serially with `--cores 1`
- avoid running bootstrap and jackknife together with other heavy workflow blocks
- the current implementation uses memory-mapped `float32` unit matrices to reduce RAM pressure

### 13. High-k Sketch Module

Input:

- `unmasked` or `masked` genomes

Actions:

- compute MinHash-style signatures for larger `k`
- estimate pairwise distances via signature Jaccard
- infer trees from sketch distance matrices

Main outputs:

- `results/sketch/signatures/{dataset}/k{k}/{accession}.tsv`
- `results/sketch/distances/{dataset}/k{k}/minhash_jaccard.tsv`
- `results/sketch/trees/{dataset}/k{k}/{method}.nwk`

## Configuration

Main configuration file:

- `config/config.yaml`

Current configurable sections:

- `metadata`
- `downloads`
- `preprocessing`
- `organelle_screen`
- `repeat_annotation`
- `kmers`
- `distances`
- `trees`
- `resampling`
- `sketch`

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

To run the broader downstream pipeline after the pre-kmer milestone:

```bash
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 4 full_analysis
```

## Current Status

Implemented:

- accession-driven metadata resolution
- genome download
- QC
- preprocessing
- organelle screening and filtering
- low-complexity annotation and masking
- pre-kmer summary reporting
- small-k spectra
- distance matrices
- tree inference
- tree comparison reports
- bootstrap and jackknife support summaries
- high-k sketch distances and trees

Not implemented yet:

- tree visualization figures
- full TE-aware repeat annotation backend

## Notes

- Samples are keyed by assembly accession.
- Metadata are resolved from NCBI E-utilities with direct accession lookup.
- Downloaded genome filenames are normalized as `data/genomes/{accession}.fna.gz`.
- The current pre-kmer workflow filters organellar contigs but does not perform general decontamination for symbionts or other non-target contaminants.
- At this stage the workflow assumes the deposited nuclear assemblies are otherwise biologically clean enough for downstream comparative analyses.
- In downstream k-mer analyses, the `unmasked` dataset refers to organelle-filtered genomes, not raw preprocessed assemblies.
- The resampling module should be launched conservatively, ideally with `--cores 1`, because it still processes large unit matrices even though they are memory-mapped.
- The current machine has a broken `/usr/local/bin/snakemake`; use `/home/bosi/miniforge3/envs/ampwrap/bin/snakemake`.

## Repository Scope

Track code, configuration and lightweight examples in Git.

Do not track:

- workflow outputs in `results/`
- Snakemake runtime state in `.snakemake/`
- downloaded genomes in `data/genomes/`
- large local resources in `resources/`
