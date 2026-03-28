import csv
from pathlib import Path

from Bio import Phylo
from phylo_utils import infer_tree_from_distance_matrix


input_path = Path(snakemake.input[0])
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
method = snakemake.params.method

with input_path.open() as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    names = header[1:]
    full_matrix = []
    row_names = []
    for row in reader:
        row_names.append(row[0])
        full_matrix.append([float(value) for value in row[1:]])

if names != row_names:
    raise ValueError("Distance matrix row and column labels do not match")

tree = infer_tree_from_distance_matrix(names, full_matrix, method)

Phylo.write(tree, str(output_path), "newick")
