# Project Roadmap

This file tracks the current implementation status of the workflow and the next milestones.

## Completed

- [x] Initialize Git repository and connect remote
- [x] Define the accession-driven project structure
- [x] Resolve assembly metadata from NCBI
- [x] Download genome FASTA files from NCBI
- [x] Implement initial assembly QC
- [x] Implement genome preprocessing
- [x] Implement organelle screening with similarity-based detection
- [x] Remove `organelle_confident` contigs
- [x] Implement low-complexity annotation with `dustmasker`
- [x] Implement low-complexity masking
- [x] Implement pre-kmer summary and Markdown report
- [x] Add smoke tests for the pre-kmer milestone
- [x] Validate the pre-kmer milestone on real genomes
- [x] Commit the pre-kmer milestone

## Small-k Phylogenomics

- [x] Confirm the correct `unmasked` dataset input for k-mers
- [x] Stabilize `workflow/scripts/compute_kmer_spectrum.py`
- [x] Generate spectra for multiple `k` values
- [x] Build feature matrices per dataset and per `k`
- [x] Validate cosine distance
- [x] Validate Jensen-Shannon distance
- [x] Produce distance matrices for `unmasked` and `masked`
- [x] Validate `nj`
- [x] Validate `upgma`
- [x] Produce Newick outputs for all supported parameter combinations
- [x] Regenerate `tree_manifest.tsv`
- [x] Regenerate `tree_comparisons.tsv`
- [x] Validate `masked` vs `unmasked` comparisons
- [x] Run the workflow on real data through k-mers, distances, trees and comparisons
- [x] Commit the small-k milestone

## Resampling

- [x] Integrate bootstrap resampling on windows
- [x] Integrate jackknife resampling on contigs
- [x] Add smoke tests for resampling
- [x] Refactor resampling to low-memory binary unit matrices
- [x] Validate minimal real-data bootstrap
- [x] Validate minimal real-data jackknife
- [x] Add conservative runtime guidance and memory guardrails
- [x] Commit the resampling milestone

## High-k Sketch

- [x] Integrate the MinHash sketch module
- [x] Compute sketch-based distance matrices
- [x] Infer sketch-based trees
- [x] Add smoke tests for sketch mode
- [x] Validate sketch mode on a minimal real-data case
- [ ] Validate the full sketch grid on real data
- [ ] Commit the sketch milestone

## Documentation Cleanup

- [ ] Ensure README fully reflects the current milestone structure
- [ ] Document clearly that `unmasked` means `organelle-filtered`
- [ ] Document explicitly that general decontamination of symbionts is not yet implemented
- [ ] Document that resampling should be run conservatively, ideally with `--cores 1`

## Future Biology-Oriented Improvements

- [ ] Affiancare the current `dustmasker` low-complexity masking with a biologically stronger repeat workflow, ideally `RepeatModeler/RepeatMasker` or an equivalent TE-aware backend

## Final Validation

- [ ] Run the full workflow on a small real dataset including sketch mode
- [ ] Check that outputs are biologically sensible across modules
- [ ] Decide whether additional contamination screening should become a formal milestone
