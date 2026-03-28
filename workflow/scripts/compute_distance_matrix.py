import csv
from pathlib import Path

from phylo_utils import compute_full_distance_matrix


metric = snakemake.params.metric
input_path = Path(snakemake.input[0])
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with input_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    feature_names = [name for name in reader.fieldnames if name != "accession"]
    data = [(row["accession"], [float(row[name]) for name in feature_names]) for row in reader]

accessions = [item[0] for item in data]
vectors = [item[1] for item in data]
full_matrix = compute_full_distance_matrix(vectors, metric)

with output_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession"] + accessions)
    for accession_a, row_values in zip(accessions, full_matrix):
        row = [accession_a]
        for distance in row_values:
            row.append(f"{distance:.12f}")
        writer.writerow(row)
