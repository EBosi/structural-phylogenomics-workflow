import csv
from pathlib import Path

from io_utils import read_fasta


def fasta_stats(path):
    n_sequences = 0
    total_bases = 0
    for _, sequence in read_fasta(path):
        n_sequences += 1
        total_bases += len(sequence)
    return n_sequences, total_bases


accession = snakemake.params.accession
raw_sequences, raw_bases = fasta_stats(snakemake.input.raw)
processed_sequences, processed_bases = fasta_stats(snakemake.input.processed)

row = {
    "sample": accession,
    "raw_sequences": raw_sequences,
    "processed_sequences": processed_sequences,
    "retained_sequences": processed_sequences,
    "removed_sequences": raw_sequences - processed_sequences,
    "retained_sequence_fraction": round((processed_sequences / raw_sequences), 6) if raw_sequences else 0,
    "raw_bases": raw_bases,
    "processed_bases": processed_bases,
    "retained_bases": processed_bases,
    "removed_bases": raw_bases - processed_bases,
    "retained_base_fraction": round((processed_bases / raw_bases), 6) if raw_bases else 0,
}

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
