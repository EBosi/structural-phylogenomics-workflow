import csv
from pathlib import Path

from io_utils import read_fasta


def compute_gc(sequence):
    gc = sum(1 for base in sequence if base in {"G", "C", "g", "c"})
    atgc = sum(1 for base in sequence if base.upper() in {"A", "C", "G", "T"})
    return (gc / atgc) * 100 if atgc else 0.0


def n50(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            return length
    return 0


sample = snakemake.params.sample
fasta_path = snakemake.input[0]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

lengths = []
total_bases = 0
total_n = 0
gc_weighted_numerator = 0.0

for _, sequence in read_fasta(fasta_path):
    seq_len = len(sequence)
    lengths.append(seq_len)
    total_bases += seq_len
    total_n += sum(1 for base in sequence if base.upper() == "N")
    gc_weighted_numerator += compute_gc(sequence) * seq_len

row = {
    "sample": sample,
    "input_fasta": fasta_path,
    "n_sequences": len(lengths),
    "total_bases": total_bases,
    "min_length": min(lengths) if lengths else 0,
    "max_length": max(lengths) if lengths else 0,
    "mean_length": round(total_bases / len(lengths), 2) if lengths else 0,
    "n50": n50(lengths),
    "gc_percent": round(gc_weighted_numerator / total_bases, 4) if total_bases else 0,
    "n_fraction_percent": round((total_n / total_bases) * 100, 4) if total_bases else 0,
}

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
