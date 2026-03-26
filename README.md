# K-mer Phylogeny Workflow

Snakemake workflow scaffold for the first preprocessing stages of an alignment-free phylogenomics pipeline on eukaryotic genomes.

Current modules:

- dataset table representation
- QC
- genome preprocessing
- repeat annotation
- repeat masking

The default repeat backend uses `dustmasker` to make the workflow runnable without `RepeatMasker`. This is useful as an MVP for low-complexity/repeat-sensitive preprocessing, but it is not a substitute for a curated TE annotation workflow.

## Repository scope

Track code, configuration and lightweight examples in Git.

Do not track:

- workflow outputs in `results/`
- Snakemake runtime state in `.snakemake/`
- large private input genomes in `data/genomes/`
- custom databases or large resources in `resources/`

## Layout

- `data/genomes/`: input assemblies (`.fa`, `.fasta`, `.fna`, optionally gzipped)
- `config/config.yaml`: workflow configuration
- `workflow/rules/`: Snakemake rules split by module
- `workflow/scripts/`: Python helper scripts
- `results/`: generated outputs

## Outputs

- `results/dataset/genome_dataset.tsv`
- `results/qc/{sample}.tsv`
- `results/qc/qc_summary.tsv`
- `results/preprocessed/{sample}.fa`
- `results/repeats/annotation/{sample}.intervals.txt`
- `results/repeats/annotation/{sample}.summary.tsv`
- `results/repeats/repeat_annotation_summary.tsv`
- `results/repeats/masked/{sample}.fa`

## Run

```bash
cd /home/bosi/kmer_phylo_workflow
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake -n
/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 4
```

## Notes

- Sample names are derived from FASTA filenames.
- Preprocessing removes short contigs and optionally excludes organelle-like sequences by header keyword.
- QC is computed on raw genomes before preprocessing.
- Repeat annotation and masking are run on preprocessed genomes.
- To switch later to `RepeatModeler/RepeatMasker`, replace the `repeats.smk` shell commands and keep the same output contract.
- In the current machine, `/usr/local/bin/snakemake` is misconfigured; use the binary in `miniforge3/envs/ampwrap/bin/`.
