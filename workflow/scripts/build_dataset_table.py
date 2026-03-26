import csv
from pathlib import Path


def detect_compression(path_str):
    return "gz" if path_str.endswith(".gz") else "plain"


rows = []
for sample, fasta in sorted(snakemake.params.genomes.items()):
    path = Path(fasta)
    rows.append(
        {
            "sample": sample,
            "assembly_path": str(path),
            "filename": path.name,
            "compression": detect_compression(str(path)),
        }
    )

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["sample", "assembly_path", "filename", "compression"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
