import csv
from pathlib import Path

from io_utils import read_fasta


def load_sequence_lengths(path):
    return {header.split()[0]: len(sequence) for header, sequence in read_fasta(path)}


def parse_interval_file(path):
    current_seq = None
    masked_bases = 0
    interval_count = 0

    with Path(path).open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq = line[1:].split()[0]
                continue
            if current_seq is None:
                continue
            if "-" not in line:
                continue
            start_str, end_str = [part.strip() for part in line.split("-", maxsplit=1)]
            start = int(start_str)
            end = int(end_str)
            interval_count += 1
            masked_bases += max(0, end - start + 1)

    return interval_count, masked_bases


sample = snakemake.params.sample
seq_lengths = load_sequence_lengths(snakemake.input.fasta)
interval_count, masked_bases = parse_interval_file(snakemake.input.intervals)
total_bases = sum(seq_lengths.values())

row = {
    "sample": sample,
    "n_sequences": len(seq_lengths),
    "total_bases": total_bases,
    "masked_intervals": interval_count,
    "masked_bases": masked_bases,
    "masked_fraction_percent": round((masked_bases / total_bases) * 100, 4) if total_bases else 0,
}

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
