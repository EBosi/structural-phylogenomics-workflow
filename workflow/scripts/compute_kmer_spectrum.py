import csv
from collections import Counter
from pathlib import Path

from io_utils import read_fasta


VALID_BASES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical_kmer(kmer):
    rc = reverse_complement(kmer)
    return min(kmer, rc)


input_fasta = snakemake.input[0]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

accession = snakemake.params.accession
dataset = snakemake.params.dataset
k = int(snakemake.params.k)
use_canonical = bool(snakemake.params.canonical)

counts = Counter()
total_valid = 0

for _, sequence in read_fasta(input_fasta):
    sequence = sequence.upper()
    if len(sequence) < k:
        continue
    for start in range(len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if set(kmer) - VALID_BASES:
            continue
        if use_canonical:
            kmer = canonical_kmer(kmer)
        counts[kmer] += 1
        total_valid += 1

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["accession", "dataset", "k", "kmer", "count", "frequency"],
        delimiter="\t",
    )
    writer.writeheader()
    for kmer in sorted(counts):
        count = counts[kmer]
        writer.writerow(
            {
                "accession": accession,
                "dataset": dataset,
                "k": k,
                "kmer": kmer,
                "count": count,
                "frequency": f"{(count / total_valid) if total_valid else 0.0:.12f}",
            }
        )
