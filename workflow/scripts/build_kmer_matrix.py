import csv
from collections import defaultdict
from pathlib import Path


input_paths = [Path(path) for path in snakemake.input]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

value_column = snakemake.params.value_column
requested_accessions = list(snakemake.params.accessions)

rows_by_accession = defaultdict(dict)
all_kmers = set()

for path in input_paths:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row["accession"]
            kmer = row["kmer"]
            rows_by_accession[accession][kmer] = row[value_column]
            all_kmers.add(kmer)

ordered_kmers = sorted(all_kmers)
fieldnames = ["accession"] + ordered_kmers

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for accession in requested_accessions:
        row = {"accession": accession}
        spectrum = rows_by_accession.get(accession, {})
        for kmer in ordered_kmers:
            row[kmer] = spectrum.get(kmer, "0")
        writer.writerow(row)
