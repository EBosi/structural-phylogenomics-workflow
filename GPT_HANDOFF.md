# GPT Handoff

## Repo

- Path: `/home/bosi/kmer_phylo_workflow`
- Current branch: `reads-mode`
- `main` and `origin/main` currently point to commit `87f73c2` (`Add configurable repeat masking backends`)

## Current Git State

The `reads-mode` branch contains uncommitted work.

Modified tracked files:

- `README.md`
- `Snakefile`
- `config/config.yaml`
- `tests/test_workflow_smoke.py`
- `workflow/rules/pre_kmer_reports.smk`
- `workflow/rules/repeats.smk`
- `workflow/scripts/build_pre_kmer_summary.py`
- `workflow/scripts/compute_sketch_distance_matrix.py`
- `workflow/scripts/io_utils.py`
- `workflow/scripts/render_pre_kmer_report.py`
- `workflow/scripts/repeat_utils.py`
- `workflow/scripts/run_repeat_annotation.py`
- `workflow/scripts/summarize_repeat_annotation.py`

New untracked files:

- `GPT_HANDOFF.md`
- `metadata/local_reads.tsv`
- `metadata/repeatmasker.tsv`
- `workflow/rules/reads.smk`
- `workflow/scripts/build_reads_manifest.py`
- `workflow/scripts/compute_filtered_read_signature.py`
- `workflow/scripts/compute_read_kmer_histogram.py`
- `workflow/scripts/compute_read_qc.py`
- `workflow/scripts/infer_read_abundance_band.py`
- `workflow/scripts/infer_tree_safe.py`
- `workflow/scripts/model_read_k_response.py`
- `workflow/scripts/read_histogram_model.py`
- `workflow/scripts/render_reads_report.py`
- `workflow/scripts/run_repeatmodeler.py`
- `workflow/scripts/score_read_samples.py`
- `workflow/scripts/select_dataset_k.py`
- `workflow/scripts/subset_fasta_prefix.py`
- `workflow/scripts/summarize_read_histograms.py`
- `workflow/scripts/summarize_read_models.py`
- `workflow/scripts/summarize_read_qc.py`

Generated untracked directories from manual/failed RepeatMasker probes:

- `RM_2.WedApr290802382026/`
- `RM_42.WedApr290808252026/`

Do not delete them silently unless the user approves cleanup; they are not part of the workflow contract.

Nothing in the reads branch has been committed yet.

## What Exists on `main`

The assembly-based workflow is mature and already committed/pushed.

Milestones already on `main`:

- `b2c1005` Bootstrap scaffold
- `158467f` Accession-driven metadata/download refactor
- `eaf0905` Pre-kmer assembly milestone
- `b0de40f` Small-k workflow
- `fd6c6fc` Low-memory bootstrap/jackknife
- `5f2a62a` High-k sketch workflow
- `aeffb62` Local genome input support
- `48f9404` Documentation cleanup
- `87f73c2` Configurable repeat masking backends

## Assembly Workflow Summary

### Inputs

- `metadata/accessions.txt`
- `metadata/local_genomes.tsv`

### Assembly pipeline

- NCBI metadata resolution
- genome download
- raw QC
- preprocessing
- organelle screening/removal
- repeat/low-complexity masking
- small-k spectra
- distance matrices
- tree inference
- tree comparison
- resampling
- high-k sketch

### Current repeat backends

- `dustmasker`
- `repeatmasker`
- `dustmasker+repeatmasker`
- `repeatmodeler` as a configured workflow component, but not yet runnable in the current env because `RepeatModeler`/`BuildDatabase` are not installed

Config location:

- `config/config.yaml`

RepeatMasker notes:

- RepeatMasker is installed in env `ampwrap`
- RepeatMasker version observed: `4.2.2`
- RMBLAST engine observed: `2.14.1+`
- Dfam version observed: `3.9`
- Dfam partitions `14` and `15` were downloaded manually into:
  - `/home/bosi/miniforge3/envs/ampwrap/share/RepeatMasker/Libraries/famdb`
- `famdb.py` shebang was patched outside the repo to use env python:
  - `/home/bosi/miniforge3/envs/ampwrap/share/RepeatMasker/famdb.py`

## Repeat Layer Roadmap / Current Status

### Point 1: Dustmasker baseline

Status: complete.

The default assembly workflow remains `repeat_annotation.backend: "dustmasker"`.

Validated baseline masked fractions:

- `GCA_000934665.2`: `3.1922%`
- `GCA_021134715.1`: `5.2863%`
- `GCF_943734735.2`: `4.4611%`
- `GCA_052721305.1`: `3.2429%`
- `protura_purged`: `8.2137%`

This is technically useful as low-complexity filtering but biologically weak for eukaryotic TE masking.

### Point 2: RepeatMasker known-library pilot

Status: implemented as a bounded, reproducible pilot target.

New configuration:

- `repeat_annotation.repeatmasker_sample_config: "metadata/repeatmasker.tsv"`
- `repeat_annotation.repeatmasker_pilot_accessions`
- `repeat_annotation.repeatmasker_pilot_backend`
- `repeat_annotation.repeatmasker_pilot_max_bases`

New target:

- `repeatmasker_pilot`

New pilot outputs:

- `results/repeats/pilot/{accession}/input.first{max_bases}.fa`
- `results/repeats/pilot/{accession}/first{max_bases}/annotation.intervals.txt`
- `results/repeats/pilot/{accession}/first{max_bases}/annotation.details.tsv`
- `results/repeats/pilot/{accession}/first{max_bases}/summary.tsv`
- `results/repeats/pilot/{accession}/first{max_bases}/classes.tsv`
- `results/repeats/pilot/summary.tsv`
- `results/repeats/pilot/class_summary.tsv`

Why bounded pilot exists:

- A full-genome Anopheles RepeatMasker run at 1 core started as `4598` batches and was operationally too slow for an interactive validation step.
- The bounded target validates the workflow contract, parsing, interval merge, source-aware reporting, and biological class breakdown without overwriting dust-only production outputs.

Anopheles sample-specific RepeatMasker config:

- File: `metadata/repeatmasker.tsv`
- Accession: `GCF_943734735.2`
- Species: `Anopheles gambiae`
- Extra args: `-uncurated -qq`

Reasoning:

- `famdb.py` showed `Anopheles gambiae` has `16` curated ancestor families and `0` curated lineage-specific families.
- With `-uncurated`, Dfam exposes `2537` lineage-specific Anopheles families.
- Curated-only masking is therefore too weak for this pilot.

Important bug fixed:

- RepeatMasker can fail on long FASTA identifiers.
- `workflow/scripts/run_repeat_annotation.py` now writes RepeatMasker-safe temporary FASTA IDs and maps parsed `.out` records back to original sequence IDs.

Pilot result on the first `1,000,000` bp of `GCF_943734735.2`:

- combined masked fraction: `21.5265%`
- dustmasker: `8.596%`
- RepeatMasker known: `18.6813%`
- RepeatMasker custom/de novo: `0%`

Observed class breakdown includes DNA transposons, LINE, LTR, RC/Helitron, SINE, Satellite, Simple_repeat, Low_complexity, and Unknown.

This result is not a whole-genome estimate. It proves the TE-aware path is functional and interpretable.

### Point 3: De novo repeat library

Status: not complete.

`workflow/scripts/run_repeatmodeler.py` and Snakemake wiring exist, but `RepeatModeler` and `BuildDatabase` are not installed in the current `ampwrap` env. Publication-grade masking for non-model eukaryotes still needs:

- RepeatModeler installation and validation.
- Per-sample de novo libraries.
- RepeatMasker run with known Dfam libraries and de novo libraries.
- Explicit reporting of known vs de novo vs DUST contributions.

### Repeat Layer Validation Commands

Recent validation:

- `python -m py_compile workflow/scripts/*.py`
- `pytest -q` -> `12 passed`
- `/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 1 --dry-run repeatmasker_pilot`
- `/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 1 repeatmasker_pilot`

`pre_kmer` dry-run remains dust-only by default, but Snakemake reports jobs as rerunnable because the repeat scripts/params changed.

## Assembly Results Already Present

Main results root:

- `/home/bosi/kmer_phylo_workflow/results`

Useful files:

- `results/reports/pre_kmer_report.md`
- `results/reports/tree_manifest.tsv`
- `results/reports/tree_comparisons.tsv`

Reference tree discussed in the thread:

- `results/trees/masked/k7/jensen_shannon/nj.nwk`

Newick:

```text
((GCA_000934665.2:0.07472,(GCF_943734735.2:0.09880,GCA_021134715.1:0.05074)Inner1:0.01511)Inner2:0.01282,GCA_052721305.1:0.08285,protura_purged:0.06756)Inner3:0.00000;
```

Species mapping used:

- `GCA_000934665.2` -> `Catajapyx aquilonaris`
- `GCA_021134715.1` -> `Daphnia pulex`
- `GCF_943734735.2` -> `Anopheles gambiae`
- `GCA_052721305.1` -> `Folsomia candida`
- `protura_purged` -> `Protura sp.`

Small-k values used:

- `5, 6, 7, 8`

Resampling k:

- `7`

Sketch k values:

- `15, 17, 21`

## Reads-Mode Prototype: Intent

The user wants a robust and potentially innovative read-based mode, not just a trivial Mash-on-reads wrapper.

Target idea:

- adaptive abundance-filtered alignment-free phylogenomics from reads
- per-sample QC and support scoring
- dataset-level K selection
- sample flagging/exclusion if quality is poor
- informative-band filtering to suppress:
  - low-frequency errors
  - high-copy repeats/organelle-heavy signal
- filtered sketching from reads

## Reads-Mode Prototype: What Has Been Implemented Locally

### New input

- `metadata/local_reads.tsv`

Columns currently expected:

- `sample_id`
- `species`
- `read1`
- `read2`
- `estimated_genome_size_bp`
- `platform`
- `notes`

### New target

- `reads_analysis`

### New reads pipeline pieces

#### Metadata + QC

- `workflow/rules/reads.smk`
- `workflow/scripts/build_reads_manifest.py`
- `workflow/scripts/compute_read_qc.py`
- `workflow/scripts/summarize_read_qc.py`

Outputs:

- `results/reads/metadata/samples.tsv`
- `results/reads/qc/{sample}.tsv`
- `results/reads/qc/summary.tsv`

#### Histogram stage

- `workflow/scripts/compute_read_kmer_histogram.py`
- `workflow/scripts/summarize_read_histograms.py`

Behavior:

- samples a bounded number of reads per file
- computes canonical k-mer counts for candidate K values
- writes count-of-counts histogram rows, not just one summary row per K

Outputs:

- `results/reads/histograms/{sample}.tsv`
- `results/reads/histograms/summary.tsv`

Candidate K values from config:

- `15, 17, 19, 21`

#### K response / sample scoring

- `workflow/scripts/model_read_k_response.py`
- `workflow/scripts/score_read_samples.py`
- `workflow/scripts/select_dataset_k.py`

Current logic:

- coverage estimated from total read bases / user-provided genome size
- histogram-derived fields include:
  - `singleton_fraction`
  - `retained_fraction`
  - `signal_peak_abundance`
  - `low_band_suggestion`
  - `high_band_suggestion`
  - `high_copy_fraction`
- sample gets:
  - `max_supported_k`
  - `sample_quality_score`
  - `recommended_action`

Outputs:

- `results/reads/models/{sample}.tsv`
- `results/reads/models/summary.tsv`
- `results/reads/sample_scores.tsv`
- `results/reads/dataset_k_selection.tsv`

#### Band inference

- `workflow/scripts/infer_read_abundance_band.py`

Output:

- `results/reads/bands/{sample}.tsv`

Current idea:

- use the selected K
- infer a low boundary near the valley after the error peak
- infer a high boundary where the count-of-counts tail becomes sparse

This is still heuristic/peak-aware, not a formal mixture fit.

#### Filtered sketch from reads

- `workflow/scripts/compute_filtered_read_signature.py`
- reused `workflow/scripts/compute_sketch_distance_matrix.py`
- new `workflow/scripts/infer_tree_safe.py`

Outputs:

- `results/reads/sketch/signatures/{sample}.tsv`
- `results/reads/sketch/distances/minhash_jaccard.tsv`
- `results/reads/sketch/tree.nwk`

Current sketch behavior:

- uses dataset-level recommended K
- resamples reads
- keeps only k-mers with abundance inside inferred sample band
- MinHash sketch built from retained k-mers

`infer_tree_safe.py` exists because the reads target may run with 0 or 1 sample and standard tree inference is not safe in that case.

### Supporting file changes

- `Snakefile` now includes `workflow/rules/reads.smk` and target `reads_analysis`
- `config/config.yaml` now has a `reads:` block
- `workflow/scripts/io_utils.py` now also provides `read_fastq()`
- `README.md` mentions the reads-mode scaffold
- `workflow/scripts/compute_sketch_distance_matrix.py` was generalized to accept either `params.accessions` or `params.samples`

## Reads-Mode Current Config

See `config/config.yaml`.

Relevant current keys:

- `metadata.local_reads_file`
- `reads.candidate_k_values`
- `reads.preferred_k_values`
- `reads.min_supported_k`
- `reads.histogram_max_reads_per_file`
- `reads.histogram_max_count_bin`
- `reads.sketch_max_reads_per_file`
- `reads.sketch_num_hashes`
- coverage thresholds under `reads.*coverage_threshold`

## Reads-Mode Validation Already Done

Validation has been on empty metadata / scaffold behavior, not on real FASTQ data.

Confirmed:

- `python3 -m py_compile` on the new reads scripts passes
- `/home/bosi/miniforge3/envs/ampwrap/bin/snakemake --cores 1 reads_analysis` completes successfully with empty `metadata/local_reads.tsv`

Produced outputs in that empty-sample state:

- `results/reads/metadata/samples.tsv`
- `results/reads/qc/summary.tsv`
- `results/reads/histograms/summary.tsv`
- `results/reads/models/summary.tsv`
- `results/reads/sample_scores.tsv`
- `results/reads/dataset_k_selection.tsv`
- `results/reads/sketch/distances/minhash_jaccard.tsv`
- `results/reads/sketch/tree.nwk`
- `results/reads/report.md`

## Important Caveat

The assistant explicitly said the reads-mode is **not ready** and is still a prototype.

Reason:

- no real FASTQ validation yet
- no proper peak fitting
- no trimming/error model integration
- histograms are sampled
- abundance-band inference is still heuristic
- no benchmark against Mash/Skmer/kWIP yet

## Recommended Next Steps

1. Commit the current `reads-mode` prototype on its own branch as an explicit prototype milestone.
2. Test the reads pipeline on real FASTQ datasets.
3. Inspect:
   - histogram shape
   - selected K
   - inferred low/high abundance bands
   - retained k-mer counts
   - whether poor samples are flagged correctly
4. Replace heuristic band inference with a more formal peak model.
5. Compare read-based distances/trees against:
   - assembly-based workflow results
   - Mash/Skmer/kWIP-type baselines

## Suggested Prompt To Another Model

Use this repo:

- `/home/bosi/kmer_phylo_workflow`

Start from branch:

- `reads-mode`

Important:

- `main` is stable for assembly-based workflow
- `reads-mode` contains uncommitted prototype work for adaptive read-based phylogenomics
- the goal is to strengthen the read-based method, especially histogram modeling, abundance-band inference, and robust K selection
- do not revert the uncommitted reads changes
