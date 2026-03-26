import csv
from pathlib import Path


input_paths = [Path(path) for path in snakemake.input]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

rows = []
fieldnames = None

for path in input_paths:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            continue
        if fieldnames is None:
            fieldnames = reader.fieldnames
        elif fieldnames != reader.fieldnames:
            raise ValueError(f"Header mismatch while merging {path}")
        rows.extend(reader)

fieldnames = fieldnames or []

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
