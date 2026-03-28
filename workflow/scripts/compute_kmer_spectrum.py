import csv
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from io_utils import read_fasta


BITS_TO_BASE = "ACGT"
INVALID = 255
CHUNK_SIZE = 2_000_000

BYTE_TO_BITS = np.full(256, INVALID, dtype=np.uint8)
BYTE_TO_BITS[ord("A")] = 0
BYTE_TO_BITS[ord("C")] = 1
BYTE_TO_BITS[ord("G")] = 2
BYTE_TO_BITS[ord("T")] = 3


def decode_kmer(code, k):
    chars = ["A"] * k
    for idx in range(k - 1, -1, -1):
        chars[idx] = BITS_TO_BASE[code & 0b11]
        code >>= 2
    return "".join(chars)


def count_chunk(chunk_bytes, k, use_canonical, counts):
    if len(chunk_bytes) < k:
        return 0

    encoded = BYTE_TO_BITS[np.frombuffer(chunk_bytes, dtype=np.uint8)]
    windows = sliding_window_view(encoded, k)
    valid = np.all(windows < 4, axis=1)
    if not np.any(valid):
        return 0

    forward = np.zeros(windows.shape[0], dtype=np.uint32)
    for idx in range(k):
        forward = (forward << 2) | windows[:, idx].astype(np.uint32)

    if use_canonical:
        reverse = np.zeros(windows.shape[0], dtype=np.uint32)
        for idx in range(k):
            reverse = (reverse << 2) | (3 - windows[:, k - 1 - idx]).astype(np.uint32)
        codes = np.minimum(forward, reverse)[valid]
    else:
        codes = forward[valid]

    bincount = np.bincount(codes, minlength=counts.shape[0])
    counts += bincount.astype(counts.dtype, copy=False)
    return int(valid.sum())


input_fasta = snakemake.input[0]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

accession = snakemake.params.accession
dataset = snakemake.params.dataset
k = int(snakemake.params.k)
use_canonical = bool(snakemake.params.canonical)

counts = np.zeros(1 << (2 * k), dtype=np.uint64)
total_valid = 0

for _, sequence in read_fasta(input_fasta):
    sequence_bytes = sequence.upper().encode("ascii")
    max_start = len(sequence_bytes) - k + 1
    if max_start <= 0:
        continue

    for start in range(0, max_start, CHUNK_SIZE):
        end = min(len(sequence_bytes), start + CHUNK_SIZE + k - 1)
        total_valid += count_chunk(sequence_bytes[start:end], k, use_canonical, counts)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["accession", "dataset", "k", "kmer", "count", "frequency"],
        delimiter="\t",
    )
    writer.writeheader()
    for code, count in enumerate(counts):
        if count == 0:
            continue
        writer.writerow(
            {
                "accession": accession,
                "dataset": dataset,
                "k": k,
                "kmer": decode_kmer(code, k),
                "count": int(count),
                "frequency": f"{(count / total_valid) if total_valid else 0.0:.12f}",
            }
        )
