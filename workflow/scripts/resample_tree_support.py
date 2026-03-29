import csv
import math
import random
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import Phylo

from phylo_utils import compute_full_distance_matrix, infer_tree_from_distance_matrix, tree_split_set


def load_matrices(matrix_paths, meta_paths):
    metadata_by_accession = {}
    for path in meta_paths:
        with Path(path).open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            row = next(reader)
            metadata_by_accession[row["accession"]] = row

    accession_units = {}
    n_features = None
    for path in matrix_paths:
        accession = Path(path).name.removesuffix(".npy")
        if accession.endswith(".units"):
            accession = accession[:-6]
        meta = metadata_by_accession.get(accession)
        if meta is None:
            raise ValueError(f"Missing metadata for {accession}")
        matrix = np.load(path, mmap_mode="r")
        if matrix.dtype != np.float32:
            raise ValueError(f"Expected float32 matrix for {accession}, found {matrix.dtype}")
        if n_features is None:
            n_features = int(meta["n_features"])
        elif n_features != int(meta["n_features"]):
            raise ValueError("Inconsistent feature dimensions across unit matrices")
        accession_units[accession] = matrix
    return accession_units, n_features or 0


def sample_unit_weights(n_units, mode, fraction, rng):
    if mode == "bootstrap":
        sampled = [rng.randrange(n_units) for _ in range(n_units)]
        return np.bincount(sampled, minlength=n_units).astype(np.float32, copy=False)
    if mode == "jackknife":
        n = max(1, math.ceil(n_units * fraction))
        weights = np.zeros(n_units, dtype=np.float32)
        weights[rng.sample(range(n_units), n)] = 1.0
        return weights
    raise ValueError(f"Unsupported resampling mode: {mode}")


def aggregate_units(matrix, weights, n_features):
    if weights.size == 0 or float(weights.sum()) == 0.0:
        return np.zeros(n_features, dtype=np.float32)
    aggregated = np.asarray(weights @ matrix, dtype=np.float32)
    total = float(aggregated.sum())
    if total == 0.0:
        return np.zeros(n_features, dtype=np.float32)
    aggregated /= total
    return aggregated


base_tree = Phylo.read(snakemake.input.tree, "newick")
base_splits = tree_split_set(base_tree)
accession_units, n_features = load_matrices(snakemake.input.matrices, snakemake.input.metadata)

accessions = sorted(accession_units)
mode = snakemake.params.mode
metric = snakemake.params.metric
method = snakemake.params.method
replicates = int(snakemake.params.replicates)
fraction = float(snakemake.params.fraction)
rng = random.Random(int(snakemake.params.seed))

support_counter = Counter()

for _ in range(replicates):
    vectors = []
    for accession in accessions:
        unit_matrix = accession_units[accession]
        weights = sample_unit_weights(unit_matrix.shape[0], mode, fraction, rng)
        vectors.append(aggregate_units(unit_matrix, weights, n_features).tolist())
    full_matrix = compute_full_distance_matrix(vectors, metric)
    tree = infer_tree_from_distance_matrix(accessions, full_matrix, method)
    replicate_splits = tree_split_set(tree)
    for split in base_splits:
        if split in replicate_splits:
            support_counter[split] += 1

support_rows = []
for split in sorted(base_splits):
    supported = support_counter[split]
    support_rows.append(
        {
            "split_taxa": ",".join(split),
            "supported_replicates": supported,
            "total_replicates": replicates,
            "support_fraction": f"{(supported / replicates) if replicates else 0.0:.12f}",
        }
    )

summary_row = {
    "mode": mode,
    "dataset": snakemake.wildcards.dataset,
    "k": snakemake.wildcards.k,
    "metric": metric,
    "method": method,
    "n_taxa": len(accessions),
    "n_reference_splits": len(base_splits),
    "replicates": replicates,
    "mean_support": f"{(sum(support_counter[split] / replicates for split in base_splits) / len(base_splits)) if base_splits and replicates else 0.0:.12f}",
    "full_support_splits": sum(1 for split in base_splits if support_counter[split] == replicates),
}

Path(snakemake.output.support).parent.mkdir(parents=True, exist_ok=True)

with Path(snakemake.output.support).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["split_taxa", "supported_replicates", "total_replicates", "support_fraction"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(support_rows)

with Path(snakemake.output.summary).open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary_row.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(summary_row)
