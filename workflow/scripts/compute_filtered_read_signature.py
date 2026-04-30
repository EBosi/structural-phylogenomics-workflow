import csv
from collections import Counter
from pathlib import Path

import numpy as np

from io_utils import read_fastq


BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
MASK64 = np.uint64((1 << 64) - 1)
FIELDNAMES = [
    "accession",
    "dataset",
    "k",
    "hash_rank",
    "hash_value",
    "sampled_reads",
    "retained_kmers",
    "low_count",
    "high_count",
    "recommended_action",
    "band_confidence",
    "signature_status",
]


def iter_canonical_kmer_codes(sequence, k):
    mask = (1 << (2 * k)) - 1
    forward = 0
    reverse = 0
    valid = 0
    for char in sequence.upper():
        bits = BASE_TO_BITS.get(char)
        if bits is None:
            forward = 0
            reverse = 0
            valid = 0
            continue
        forward = ((forward << 2) | bits) & mask
        reverse = (reverse >> 2) | ((3 - bits) << (2 * (k - 1)))
        valid += 1
        if valid >= k:
            yield min(forward, reverse)


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
        return np.sort(values)
    idx = num_hashes - 1
    return np.sort(np.partition(values, idx)[:num_hashes])


def write_metadata_only(
    output_path: Path,
    sample_id: str,
    k: int,
    low_count: int,
    high_count: int,
    recommended_action: str,
    band_confidence: str,
    signature_status: str,
    sampled_reads: int = 0,
    retained_kmers: int = 0,
) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "accession": sample_id,
                "dataset": "reads_filtered",
                "k": k,
                "hash_rank": 0,
                "hash_value": "",
                "sampled_reads": sampled_reads,
                "retained_kmers": retained_kmers,
                "low_count": low_count,
                "high_count": high_count,
                "recommended_action": recommended_action,
                "band_confidence": band_confidence,
                "signature_status": signature_status,
            }
        )


manifest_path = Path(snakemake.input.manifest)
band_path = Path(snakemake.input.band)
dataset_k_path = Path(snakemake.input.dataset_k)
sample_id = snakemake.wildcards.sample
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

num_hashes = int(snakemake.params.num_hashes)
max_reads_per_file = int(snakemake.params.max_reads_per_file)
exclude_low_confidence_bands = bool(snakemake.params.exclude_low_confidence_bands)

manifest_row = None
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["sample_id"] == sample_id:
            manifest_row = row
            break
if manifest_row is None:
    raise ValueError(f"Sample {sample_id} not found in manifest {manifest_path}")

band_row = next(csv.DictReader(band_path.open(), delimiter="\t"))
dataset_row = next(csv.DictReader(dataset_k_path.open(), delimiter="\t"))

k = int(dataset_row["recommended_k"])
selected_k = int(band_row["selected_k"])
low_count = int(band_row["low_count"])
high_count = int(band_row["high_count"])
recommended_action = band_row["recommended_action"]
band_confidence = band_row["band_confidence"]
sample_supports_dataset_k = str(band_row["sample_supports_dataset_k"]) == "1"

if selected_k != k:
    raise ValueError(f"Band file selected_k={selected_k} does not match dataset recommended_k={k}")

skip_status = None
if recommended_action == "exclude":
    skip_status = "excluded_by_qc"
elif not sample_supports_dataset_k:
    skip_status = "unsupported_dataset_k"
elif exclude_low_confidence_bands and band_confidence == "low":
    skip_status = "low_confidence_band"

if skip_status is not None:
    write_metadata_only(
        output_path,
        sample_id,
        k,
        low_count,
        high_count,
        recommended_action,
        band_confidence,
        skip_status,
    )
else:
    read_paths = [Path(manifest_row["read1"])]
    if manifest_row["read2"]:
        read_paths.append(Path(manifest_row["read2"]))

    counts = Counter()
    sampled_reads = 0
    for read_path in read_paths:
        reads_seen = 0
        for _, sequence, _ in read_fastq(read_path):
            reads_seen += 1
            sampled_reads += 1
            for code in iter_canonical_kmer_codes(sequence, k):
                counts[code] += 1
            if reads_seen >= max_reads_per_file:
                break

    retained_codes = np.array(
        [code for code, abundance in counts.items() if low_count <= abundance <= high_count],
        dtype=np.uint64,
    )
    retained_hashes = retain_smallest(splitmix64_array(retained_codes), num_hashes)

    if retained_hashes.size == 0:
        write_metadata_only(
            output_path,
            sample_id,
            k,
            low_count,
            high_count,
            recommended_action,
            band_confidence,
            "no_retained_kmers",
            sampled_reads=sampled_reads,
            retained_kmers=int(retained_codes.size),
        )
    else:
        with output_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
            writer.writeheader()
            for idx, value in enumerate(retained_hashes.tolist(), start=1):
                writer.writerow(
                    {
                        "accession": sample_id,
                        "dataset": "reads_filtered",
                        "k": k,
                        "hash_rank": idx,
                        "hash_value": int(value),
                        "sampled_reads": sampled_reads,
                        "retained_kmers": int(retained_codes.size),
                        "low_count": low_count,
                        "high_count": high_count,
                        "recommended_action": recommended_action,
                        "band_confidence": band_confidence,
                        "signature_status": "ok",
                    }
                )
