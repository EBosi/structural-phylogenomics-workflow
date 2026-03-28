import csv
import subprocess
from pathlib import Path


def fetch_first_fasta(term, esearch_path, efetch_path):
    cmd = (
        f"{esearch_path} -db nucleotide -query '{term}' | "
        f"{efetch_path} -format fasta"
    )
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        text=True,
        capture_output=True,
    )
    fasta_text = result.stdout.strip()
    if not fasta_text:
        return None
    lines = fasta_text.splitlines()
    first_record = []
    seen_header = False
    for line in lines:
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
        if seen_header:
            first_record.append(line)
    return "\n".join(first_record).strip() if first_record else None


assemblies = []
with Path(snakemake.input[0]).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    assemblies.extend(reader)

esearch_path = snakemake.params.esearch
efetch_path = snakemake.params.efetch
organelle_types = list(snakemake.params.organelle_types)

seen_organisms = set()
table_rows = []
fasta_records = []

for row in assemblies:
    organism = row["organism_name"].split(" (")[0]
    if organism in seen_organisms:
        continue
    seen_organisms.add(organism)
    for organelle_type in organelle_types:
        term = f"\"{organism}\"[Organism] AND {organelle_type} AND complete genome"
        fasta_text = fetch_first_fasta(term, esearch_path, efetch_path)
        if not fasta_text:
            continue
        header = fasta_text.splitlines()[0][1:]
        table_rows.append(
            {
                "organism_name": organism,
                "organelle_type": organelle_type,
                "nuccore_id": header.split()[0],
                "fasta_header": header,
            }
        )
        fasta_records.append(fasta_text)

fasta_path = Path(snakemake.output.fasta)
table_path = Path(snakemake.output.table)
fasta_path.parent.mkdir(parents=True, exist_ok=True)
table_path.parent.mkdir(parents=True, exist_ok=True)
fasta_path.write_text("\n".join(fasta_records) + ("\n" if fasta_records else ""))

with table_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["organism_name", "organelle_type", "nuccore_id", "fasta_header"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(table_rows)
