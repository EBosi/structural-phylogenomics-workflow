#!/usr/bin/env python3
"""
Compute an exact pairwise k-mer set distance panel from FASTA inputs.

This script is intended for:
- toy validation
- small subsets
- methodological checks

It is not designed to scale to large eukaryotic assembly panels. It holds exact
unique-kmer sets in memory and should be treated as a validation utility rather
than a production workflow component.
"""

import argparse
import csv
import math
from pathlib import Path

from io_utils import read_fasta


OUTPUT_COLUMNS = [
    "sample_a",
    "sample_b",
    "n_a",
    "n_b",
    "intersection",
    "union",
    "jaccard",
    "jaccard_distance",
    "dice",
    "dice_distance",
    "overlap_coefficient",
    "overlap_distance",
    "containment_a_in_b",
    "containment_b_in_a",
    "max_containment",
    "binary_cosine",
    "binary_cosine_distance",
    "mash_distance",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", nargs="+", help="Input FASTA or FASTA.gz files.")
    parser.add_argument("--sample-name", dest="sample_names", action="append", default=[], help="Optional sample name for each FASTA, in the same order as inputs.")
    parser.add_argument("--k", type=int, required=True, help="k-mer size")
    parser.add_argument(
        "--canonical",
        choices=["true", "false"],
        default="true",
        help="Whether to canonicalize k-mers with reverse complements.",
    )
    parser.add_argument("--output", required=True, help="Output TSV path.")
    args = parser.parse_args()
    if args.k <= 0:
        parser.error("--k must be positive")
    if args.sample_names and len(args.sample_names) != len(args.fasta):
        parser.error("--sample-name must be provided exactly once per FASTA when used")
    return args


def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def canonical_kmer(kmer):
    rc = reverse_complement(kmer)
    return min(kmer, rc)


def infer_sample_name(path_str):
    name = Path(path_str).name
    for suffix in (".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path_str).stem


def iter_valid_kmers(sequence, k, canonical):
    seq = sequence.upper()
    if len(seq) < k:
        return
    for idx in range(len(seq) - k + 1):
        kmer = seq[idx : idx + k]
        if any(base not in "ACGT" for base in kmer):
            continue
        yield canonical_kmer(kmer) if canonical else kmer


def unique_kmer_set_from_fasta(path, k, canonical):
    kmers = set()
    for _, sequence in read_fasta(path):
        kmers.update(iter_valid_kmers(sequence, k, canonical))
    return kmers


def safe_ratio(numerator, denominator):
    if denominator == 0:
        return None
    return numerator / denominator


def mash_distance_from_jaccard(jaccard, k):
    if jaccard is None:
        return None
    if jaccard == 0:
        return math.inf
    return -(1.0 / k) * math.log((2.0 * jaccard) / (1.0 + jaccard))


def compute_distance_panel(set_a, set_b, k):
    n_a = len(set_a)
    n_b = len(set_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)

    jaccard = safe_ratio(intersection, union)
    dice = safe_ratio(2 * intersection, n_a + n_b)
    overlap = safe_ratio(intersection, min(n_a, n_b))
    containment_a_in_b = safe_ratio(intersection, n_a)
    containment_b_in_a = safe_ratio(intersection, n_b)
    max_containment = safe_ratio(intersection, max(n_a, n_b))
    binary_cosine = safe_ratio(intersection, math.sqrt(n_a * n_b))

    return {
        "n_a": n_a,
        "n_b": n_b,
        "intersection": intersection,
        "union": union,
        "jaccard": jaccard,
        "jaccard_distance": None if jaccard is None else 1.0 - jaccard,
        "dice": dice,
        "dice_distance": None if dice is None else 1.0 - dice,
        "overlap_coefficient": overlap,
        "overlap_distance": None if overlap is None else 1.0 - overlap,
        "containment_a_in_b": containment_a_in_b,
        "containment_b_in_a": containment_b_in_a,
        "max_containment": max_containment,
        "binary_cosine": binary_cosine,
        "binary_cosine_distance": None if binary_cosine is None else 1.0 - binary_cosine,
        "mash_distance": mash_distance_from_jaccard(jaccard, k),
    }


def format_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        if value == 0.0:
            return "0.0"
        return repr(value)
    return str(value)


def main():
    args = parse_args()
    canonical = args.canonical == "true"
    sample_names = args.sample_names or [infer_sample_name(path) for path in args.fasta]
    if len(sample_names) != len(set(sample_names)):
        duplicates = sorted({name for name in sample_names if sample_names.count(name) > 1})
        raise SystemExit(f"Duplicate sample names are not allowed: {', '.join(duplicates)}")

    samples = []
    for sample_name, fasta_path in zip(sample_names, args.fasta):
        kmers = unique_kmer_set_from_fasta(fasta_path, args.k, canonical)
        samples.append((sample_name, kmers))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample_a, set_a in samples:
            for sample_b, set_b in samples:
                panel = compute_distance_panel(set_a, set_b, args.k)
                row = {"sample_a": sample_a, "sample_b": sample_b}
                row.update({key: format_value(value) for key, value in panel.items()})
                writer.writerow(row)


if __name__ == "__main__":
    main()
