import csv
from collections import Counter
from pathlib import Path

from io_utils import read_fastq


BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}


def iter_canonical_kmer_codes(sequence, k):
    mask = (1 << (2 * k)) - 1
    forward = 0
    reverse = 0
    valid = 0

    for char in sequence.upper():
        if char not in BASE_TO_BITS:
            forward = 0
            reverse = 0
            valid = 0
            continue

        bits = BASE_TO_BITS[char]
        forward = ((forward << 2) | bits) & mask
        reverse = (reverse >> 2) | ((3 - bits) << (2 * (k - 1)))
        valid += 1
        if valid >= k:
            yield min(forward, reverse)


def count_file_kmers(path, k, max_reads):
    counts = Counter()
    reads_seen = 0
    total_kmers = 0
    for _, sequence, _ in read_fastq(path):
        reads_seen += 1
        for code in iter_canonical_kmer_codes(sequence, k):
            counts[code] += 1
            total_kmers += 1
        if reads_seen >= max_reads:
            break
    return counts, reads_seen, total_kmers


manifest = Path(snakemake.input.manifest)
sample_id = snakemake.wildcards.sample
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

candidate_k_values = [int(value) for value in snakemake.params.candidate_k_values]
max_reads_per_file = int(snakemake.params.max_reads_per_file)
max_count_bin = int(snakemake.params.max_count_bin)

record = None
with manifest.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["sample_id"] == sample_id:
            record = row
            break

if record is None:
    raise ValueError(f"Sample {sample_id} not found in manifest {manifest}")

read_paths = [Path(record["read1"])]
if record["read2"]:
    read_paths.append(Path(record["read2"]))

rows = []
for k_value in candidate_k_values:
    combined_counts = Counter()
    sampled_reads = 0
    total_kmers = 0
    for read_path in read_paths:
        counts, reads_seen, file_total_kmers = count_file_kmers(read_path, k_value, max_reads_per_file)
        combined_counts.update(counts)
        sampled_reads += reads_seen
        total_kmers += file_total_kmers

    if not combined_counts:
        rows.append(
            {
                "sample_id": sample_id,
                "k": k_value,
                "sampled_reads": sampled_reads,
                "total_kmers": total_kmers,
                "distinct_kmers": 0,
                "count_bin": 0,
                "n_distinct_kmers_at_count": 0,
            }
        )
        continue

    abundance_histogram = Counter()
    for abundance in combined_counts.values():
        abundance_histogram[min(abundance, max_count_bin)] += 1

    for count_bin in sorted(abundance_histogram):
        rows.append(
            {
                "sample_id": sample_id,
                "k": k_value,
                "sampled_reads": sampled_reads,
                "total_kmers": total_kmers,
                "distinct_kmers": len(combined_counts),
                "count_bin": count_bin,
                "n_distinct_kmers_at_count": abundance_histogram[count_bin],
            }
        )

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "k",
            "sampled_reads",
            "total_kmers",
            "distinct_kmers",
            "count_bin",
            "n_distinct_kmers_at_count",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
