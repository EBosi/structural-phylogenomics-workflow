import csv
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def build_lower_triangle(names, full_matrix):
    lower = []
    for i, _ in enumerate(names):
        lower.append(full_matrix[i][: i + 1])
    return lower


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

distance_matrix = DistanceMatrix(names=names, matrix=build_lower_triangle(names, full_matrix))
constructor = DistanceTreeConstructor()

if method == "nj":
    tree = constructor.nj(distance_matrix)
elif method == "upgma":
    tree = constructor.upgma(distance_matrix)
else:
    raise ValueError(f"Unsupported tree method: {method}")

Phylo.write(tree, str(output_path), "newick")
