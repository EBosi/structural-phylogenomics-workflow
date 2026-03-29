import csv
from pathlib import Path


def jaccard_distance(set_a, set_b):
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return 1.0 - (len(set_a & set_b) / union)


signatures = {}
for path in snakemake.input:
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        current_accession = None
        values = set()
        for row in reader:
            current_accession = row["accession"]
            values.add(int(row["hash_value"]))
        if current_accession is None:
            continue
        signatures[current_accession] = values

accessions = [acc for acc in snakemake.params.accessions if acc in signatures]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession"] + accessions)
    for acc_a in accessions:
        row = [acc_a]
        for acc_b in accessions:
            row.append(f"{jaccard_distance(signatures[acc_a], signatures[acc_b]):.12f}")
        writer.writerow(row)
