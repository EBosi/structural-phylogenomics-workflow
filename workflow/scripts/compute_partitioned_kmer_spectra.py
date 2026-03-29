import csv
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from io_utils import read_fasta


INVALID = 255
BYTE_TO_BITS = np.full(256, INVALID, dtype=np.uint8)
BYTE_TO_BITS[ord("A")] = 0
BYTE_TO_BITS[ord("C")] = 1
BYTE_TO_BITS[ord("G")] = 2
BYTE_TO_BITS[ord("T")] = 3


def reverse_complement_code(code, k):
    rc = 0
    for _ in range(k):
        base = code & 0b11
        rc = (rc << 2) | (3 - base)
        code >>= 2
    return rc


def feature_codes(k, use_canonical):
    codes = []
    for code in range(1 << (2 * k)):
        if use_canonical and code > reverse_complement_code(code, k):
            continue
        codes.append(code)
    return np.array(codes, dtype=np.uint32)


def iter_units(header, sequence, unit_type, window_size, step_size, min_window_length):
    if unit_type == "contig":
        if len(sequence) >= min_window_length:
            yield f"{header.split()[0]}:contig", sequence
        return

    if unit_type != "window":
        raise ValueError(f"Unsupported unit_type: {unit_type}")

    if len(sequence) < min_window_length:
        return

    seq_name = header.split()[0]
    for start in range(0, len(sequence), step_size):
        end = min(start + window_size, len(sequence))
        window_seq = sequence[start:end]
        if len(window_seq) < min_window_length:
            continue
        yield f"{seq_name}:{start + 1}-{end}", window_seq
        if end == len(sequence):
            break


def count_unit_frequencies(unit_seq, k, use_canonical, selected_codes):
    encoded = BYTE_TO_BITS[np.frombuffer(unit_seq.upper().encode("ascii"), dtype=np.uint8)]
    if encoded.shape[0] < k:
        return None

    windows = sliding_window_view(encoded, k)
    valid = np.all(windows < 4, axis=1)
    total_valid = int(valid.sum())
    if total_valid == 0:
        return None

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

    counts = np.bincount(codes, minlength=1 << (2 * k))
    return (counts[selected_codes] / total_valid).astype(np.float32, copy=False)


input_fasta = snakemake.input[0]
matrix_output = Path(snakemake.output.matrix)
meta_output = Path(snakemake.output.meta)
matrix_output.parent.mkdir(parents=True, exist_ok=True)

accession = snakemake.params.accession
dataset = snakemake.params.dataset
k = int(snakemake.params.k)
unit_type = snakemake.params.unit_type
window_size = int(snakemake.params.window_size)
step_size = int(snakemake.params.step_size) if int(snakemake.params.step_size) > 0 else 1
min_window_length = int(snakemake.params.min_window_length)
use_canonical = bool(snakemake.params.canonical)
selected_codes = feature_codes(k, use_canonical)

unit_ids = []
rows = []
for header, sequence in read_fasta(input_fasta):
    for unit_id, unit_seq in iter_units(header, sequence, unit_type, window_size, step_size, min_window_length):
        vector = count_unit_frequencies(unit_seq, k, use_canonical, selected_codes)
        if vector is None:
            continue
        unit_ids.append(unit_id)
        rows.append(vector)

matrix = (
    np.vstack(rows).astype(np.float32, copy=False)
    if rows
    else np.zeros((0, selected_codes.shape[0]), dtype=np.float32)
)
np.save(matrix_output, matrix)

with meta_output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["accession", "dataset", "k", "unit_type", "n_units", "n_features", "dtype"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "accession": accession,
            "dataset": dataset,
            "k": k,
            "unit_type": unit_type,
            "n_units": matrix.shape[0],
            "n_features": matrix.shape[1],
            "dtype": str(matrix.dtype),
        }
    )
