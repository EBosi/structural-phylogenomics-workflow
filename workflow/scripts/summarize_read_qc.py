import csv
from pathlib import Path


manifest_path = Path(snakemake.input.manifest)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

rows = []
for sample_report in snakemake.input.sample_reports:
    with Path(sample_report).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows.extend(reader)

rows.sort(key=lambda row: row["sample_id"])

with output_path.open("w", newline="") as handle:
    if rows:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    else:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id"])
