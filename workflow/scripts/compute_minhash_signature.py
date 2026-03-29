import csv
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from io_utils import read_fasta


INVALID = 255
CHUNK_SIZE = 2_000_000
MASK64 = np.uint64((1 << 64) - 1)

BYTE_TO_BITS = np.full(256, INVALID, dtype=np.uint8)
BYTE_TO_BITS[ord("A")] = 0
BYTE_TO_BITS[ord("C")] = 1
BYTE_TO_BITS[ord("G")] = 2
BYTE_TO_BITS[ord("T")] = 3


def splitmix64_array(values):
    values = (values + np.uint64(0x9E3779B97F4A7C15)) & MASK64
    values = ((values ^ (values >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)) & MASK64
    values = ((values ^ (values >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)) & MASK64
    return (values ^ (values >> np.uint64(31))) & MASK64


def retain_smallest(values, num_hashes):
    if values.size == 0:
        return values
    values = np.unique(values)
    if values.size <= num_hashes:
        return values
    idx = num_hashes - 1
    partitioned = np.partition(values, idx)[:num_hashes]
    return np.sort(partitioned)


def canonical_hashes_for_chunk(sequence_bytes, k, use_canonical):
    encoded = BYTE_TO_BITS[np.frombuffer(sequence_bytes, dtype=np.uint8)]
    if encoded.shape[0] < k:
        return np.array([], dtype=np.uint64)

    windows = sliding_window_view(encoded, k)
    valid = np.all(windows < 4, axis=1)
    if not np.any(valid):
        return np.array([], dtype=np.uint64)

    forward = np.zeros(windows.shape[0], dtype=np.uint64)
    for idx in range(k):
        forward = (forward << np.uint64(2)) | windows[:, idx].astype(np.uint64)

    if use_canonical:
        reverse = np.zeros(windows.shape[0], dtype=np.uint64)
        for idx in range(k):
            reverse = (reverse << np.uint64(2)) | (np.uint64(3) - windows[:, k - 1 - idx].astype(np.uint64))
        codes = np.minimum(forward, reverse)[valid]
    else:
        codes = forward[valid]

    return splitmix64_array(codes)


input_fasta = snakemake.input[0]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

accession = snakemake.params.accession
dataset = snakemake.params.dataset
k = int(snakemake.params.k)
num_hashes = int(snakemake.params.num_hashes)
use_canonical = bool(snakemake.params.canonical)

selected = np.array([], dtype=np.uint64)

for _, sequence in read_fasta(input_fasta):
    sequence_bytes = sequence.upper().encode("ascii")
    max_start = len(sequence_bytes) - k + 1
    if max_start <= 0:
        continue

    for start in range(0, max_start, CHUNK_SIZE):
        end = min(len(sequence_bytes), start + CHUNK_SIZE + k - 1)
        chunk_hashes = canonical_hashes_for_chunk(sequence_bytes[start:end], k, use_canonical)
        if chunk_hashes.size == 0:
            continue
        chunk_smallest = retain_smallest(chunk_hashes, num_hashes)
        selected = retain_smallest(np.concatenate([selected, chunk_smallest]), num_hashes)

smallest = [int(value) for value in selected.tolist()]

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["accession", "dataset", "k", "hash_rank", "hash_value"],
        delimiter="\t",
    )
    writer.writeheader()
    for idx, value in enumerate(smallest, start=1):
        writer.writerow(
            {
                "accession": accession,
                "dataset": dataset,
                "k": k,
                "hash_rank": idx,
                "hash_value": value,
            }
        )
