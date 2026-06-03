# Structural Phylogenomics Workflow

Snakemake workflow for an accession-driven, alignment-free phylogenomics pipeline on eukaryotic genomes, with optional support for locally supplied assemblies. The current implementation covers metadata resolution, genome download, QC, preprocessing, organelle screening, low-complexity masking, small-k spectrum generation, distance calculation, tree inference and tree comparison.

The default repeat backend uses `dustmasker` as a lightweight MVP. This is useful to get the workflow running, but it is not a substitute for a curated `RepeatModeler/RepeatMasker` TE annotation workflow. The repeat layer can be configured to use `dustmasker`, `RepeatMasker`, or a combined `dustmasker+repeatmasker` mode.

## Goal

The project is designed to start from a minimal user input, a list of assembly accession IDs and optionally a table of local genomes, and produce standardized metadata tables plus downstream phylogenomic outputs:

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
- `metadata/local_genomes.tsv`

One accession per line, for example:

```text
assembly_001
assembly_002
assembly_003
assembly_004
```

Alternative input from the command line:

```bash
snakemake --cores 4 --config accessions="assembly_001,assembly_002"
```

Optional local genome table:

```tsv
accession	organism_name	assembly_name	assembly_level	source_db	local_path
local_sample_001	Local species A	Local_assembly_v1	scaffold	local	data/genomes/local_sample_001.fna.gz
```

The local table must contain at least:

- `accession`
- `local_path`

## Workflow Steps

### 1. Metadata Resolution

Input:

- accession list
- optional local genome table
- NCBI E-utilities

Actions:

- resolve each accession against NCBI assembly metadata
- merge user-provided local genome records
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
- skip download for samples already provided as local FASTA files
- normalize downloaded filenames as accession-based paths

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
- or annotate repeat intervals with a TE-aware backend such as `RepeatMasker`
- or merge both interval sources in a combined mode

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

- `metadata/`: user-provided accession lists and local genome tables
- `data/genomes/`: downloaded or user-supplied assemblies
- `config/`: workflow configuration
- `workflow/rules/`: Snakemake rules by module
- `workflow/scripts/`: Python helper scripts
- `resources/`: downloaded reference tables and local resources
- `results/`: generated workflow outputs

## Portable Setup

Clone the repository on the target system:

```bash
git clone https://github.com/EBosi/structural-phylogenomics-workflow.git
cd structural-phylogenomics-workflow
```

Create and activate the workflow environment:

```bash
mamba env create -f environment.yml
mamba activate structural-phylogenomics-workflow
```

If `mamba` is not available, use Conda with the same file:

```bash
conda env create -f environment.yml
conda activate structural-phylogenomics-workflow
```

Verify that the external tools configured in `config/config.yaml` are available in the active `PATH`:

```bash
command -v curl
command -v esearch
command -v efetch
command -v blastn
command -v makeblastdb
command -v dustmasker
```

Dependency scope:

- Python smoke tests require `python`, `pytest` and the Python libraries declared in `environment.yml`.
- Snakemake dry-runs require `snakemake` and the repository files, but no NCBI downloads when no samples are configured.
- Real `pre_kmer` runs require the configured command-line tools in `PATH`: `curl`, `esearch`, `efetch`, `blastn`, `makeblastdb` and `dustmasker` for the default backend.
- `full_analysis` reuses the pre-kmer dependencies and runs the downstream Python k-mer, distance, tree, resampling and sketch scripts.
- Advanced repeat masking with `RepeatMasker` or `RepeatModeler` is optional and requires separate installation/configuration if those tools are not present in the environment.

Run the smoke tests. These tests use toy data and do not require NCBI downloads:

```bash
pytest
```

Check the Snakemake graph without running jobs:

```bash
snakemake -n
```

If `metadata/local_genomes.tsv` contains local samples, every `local_path` must point to an existing FASTA file before the dry-run can complete. For an accession-only run, leave `metadata/local_genomes.tsv` with just the header row.

To run from a directory of local assemblies, keep `metadata/accessions.txt` empty and pass the directory at runtime. FASTA files with `.fa`, `.fasta`, `.fna` and `.gz` variants are staged into `data/genomes/{sample}.fna.gz`; the sample ID is the filename without the FASTA suffix and must contain only letters, numbers, dots, underscores or hyphens.

```bash
snakemake -n --config local_genomes_dir=/path/to/assemblies
snakemake --cores 4 pre_kmer --config local_genomes_dir=/path/to/assemblies
snakemake --cores 4 full_analysis --config local_genomes_dir=/path/to/assemblies
```

Run the default pre-kmer milestone:

```bash
snakemake --cores 4 pre_kmer
```

Run the broader downstream analysis only after the pre-kmer milestone is configured and working:

```bash
snakemake --cores 4 full_analysis
```

`RepeatMasker` and `RepeatModeler` are not required for the default `dustmasker` backend. The `dustmasker` mode is a lightweight MVP for low-complexity masking. Treat `RepeatMasker` as an optional backend; if you set `repeat_annotation.backend` to `repeatmasker` or `dustmasker+repeatmasker`, install and configure `RepeatMasker` separately if it is not available in the workflow environment.

## Continuous Integration

The GitHub Actions workflow in `.github/workflows/ci.yml` is intentionally minimal. It creates the Conda/Mamba environment, runs `pytest`, and runs `snakemake -n` with no configured samples. The CI does not contact NCBI, download genomes, run BLAST jobs, run RepeatMasker, or execute `full_analysis`.

## Run

```bash
snakemake -n
snakemake --cores 4
```

To run the broader downstream pipeline after the pre-kmer milestone:

```bash
snakemake --cores 4 full_analysis
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
- automated library construction for TE-aware repeat annotation

## Notes

- Samples are keyed by a stable sample identifier in the `accession` column. For NCBI-backed samples this is the assembly accession; for local genomes it is the user-provided local sample ID.
- Metadata are resolved from NCBI E-utilities with direct accession lookup.
- Downloaded genome filenames are normalized as `data/genomes/{accession}.fna.gz`. Local genomes can be added through `metadata/local_genomes.tsv`, or by passing `--config local_genomes_dir=/path/to/assemblies` to stage a directory of local FASTA files.
- The current pre-kmer workflow filters organellar contigs but does not perform general decontamination for symbionts or other non-target contaminants.
- At this stage the workflow assumes the deposited nuclear assemblies are otherwise biologically clean enough for downstream comparative analyses.
- In downstream k-mer analyses, the `unmasked` dataset refers to organelle-filtered genomes, not raw preprocessed assemblies.
- The `repeat_annotation.backend` setting controls whether masking is driven by `dustmasker`, `RepeatMasker`, or both.
- The resampling module should be launched conservatively, ideally with `--cores 1`, because it still processes large unit matrices even though they are memory-mapped.
- Use the Snakemake executable from the intended workflow environment rather than assuming a system-wide `snakemake` binary is correct.

## Repository Scope

Track code, configuration and lightweight examples in Git.

Do not track:

- workflow outputs in `results/`
- Snakemake runtime state in `.snakemake/`
- downloaded genomes in `data/genomes/`
- large local resources in `resources/`
