# Project Roadmap

This file tracks the implementation status of the workflow and the next milestones.

## Completed

- [x] Initialize Git repository and connect remote
- [x] Define accession-driven project structure
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

## Pending Small Cleanup

- [ ] Commit the README note stating that the workflow does not yet perform general decontamination of symbionts or other non-target contaminants

## K-mer Core

- [ ] Stabilize `workflow/scripts/compute_kmer_spectrum.py`
- [x] Confirm the correct `unmasked` dataset input for k-mers
- [ ] Generate spectra for multiple `k` values
- [ ] Build feature matrices per dataset and per `k`
- [ ] Verify spectrum normalization choices
- [ ] Add smoke tests dedicated to the k-mer modules

## Distance Matrices

- [ ] Stabilize `workflow/scripts/compute_distance_matrix.py`
- [ ] Validate cosine distance
- [ ] Validate Jensen-Shannon distance
- [ ] Produce distance matrices for `unmasked` and `masked` datasets

## Tree Inference

- [ ] Stabilize `workflow/scripts/infer_tree.py`
- [ ] Validate `nj`
- [ ] Validate `upgma`
- [ ] Produce Newick outputs for all supported parameter combinations

## Tree Comparison

- [ ] Stabilize `workflow/scripts/compare_trees.py`
- [ ] Regenerate `tree_manifest.tsv`
- [ ] Regenerate `tree_comparisons.tsv`
- [ ] Validate `masked` vs `unmasked` comparisons

## Robustness

- [ ] Integrate bootstrap resampling on windows
- [ ] Integrate jackknife resampling on contigs
- [ ] Verify clade support summaries
- [ ] Add smoke tests for resampling

## High-k Sketch

- [ ] Integrate the MinHash sketch module
- [ ] Compute sketch-based distance matrices
- [ ] Infer sketch-based trees
- [ ] Add smoke tests for sketch mode

## Documentation

- [ ] Update README after the k-mer milestone
- [ ] Document clearly the current analysis datasets:
- [ ] `organelle-filtered`
- [ ] `masked`
- [ ] Document explicitly that general decontamination is not yet implemented

## Final Validation

- [ ] Run the workflow on real data through the k-mer and distance stages
- [ ] Run the full workflow on a small real dataset
- [ ] Check that outputs are biologically sensible
- [ ] Commit the k-mer milestone
