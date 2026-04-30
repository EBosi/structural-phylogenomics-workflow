import csv
from pathlib import Path

from io_utils import read_fastq

def fastq_stats(path: Path) -> dict[str, float]:
    read_count = 0
    total_bases = 0
    total_gc = 0
    total_n = 0
    total_quality = 0
    min_length = None
    max_length = 0

    for _, sequence, quality in read_fastq(path):
        read_count += 1
        sequence_upper = sequence.upper()
        seq_length = len(sequence_upper)
        total_bases += seq_length
        total_gc += sequence_upper.count("G") + sequence_upper.count("C")
        total_n += sequence_upper.count("N")
        total_quality += sum(ord(char) - 33 for char in quality)
        min_length = seq_length if min_length is None else min(min_length, seq_length)
        max_length = max(max_length, seq_length)

    mean_length = (total_bases / read_count) if read_count else 0.0
    gc_percent = (100.0 * total_gc / total_bases) if total_bases else 0.0
    n_fraction = (total_n / total_bases) if total_bases else 0.0
    mean_quality = (total_quality / total_bases) if total_bases else 0.0
    return {
        "read_count": read_count,
        "total_bases": total_bases,
        "mean_read_length": mean_length,
        "min_read_length": min_length or 0,
        "max_read_length": max_length,
        "gc_percent": gc_percent,
        "n_fraction": n_fraction,
        "mean_quality": mean_quality,
    }


manifest = Path(snakemake.input.manifest)
sample_id = snakemake.wildcards.sample
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

record = None
with manifest.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["sample_id"] == sample_id:
            record = row
            break

if record is None:
    raise ValueError(f"Sample {sample_id} not found in manifest {manifest}")

read1_path = Path(record["read1"])
read2_path = Path(record["read2"]) if record["read2"] else None
read1_stats = fastq_stats(read1_path)
read2_stats = fastq_stats(read2_path) if read2_path else None

layout = "paired-end" if read2_stats is not None else "single-end"
total_reads = read1_stats["read_count"] + (read2_stats["read_count"] if read2_stats else 0)
total_bases = read1_stats["total_bases"] + (read2_stats["total_bases"] if read2_stats else 0)
mean_read_length = total_bases / total_reads if total_reads else 0.0
pair_count_mismatch = bool(read2_stats and read1_stats["read_count"] != read2_stats["read_count"])
gc_percent = (
    (
        read1_stats["gc_percent"] * read1_stats["total_bases"]
        + ((read2_stats["gc_percent"] * read2_stats["total_bases"]) if read2_stats else 0.0)
    )
    / total_bases
    if total_bases
    else 0.0
)
n_fraction = (
    (
        read1_stats["n_fraction"] * read1_stats["total_bases"]
        + ((read2_stats["n_fraction"] * read2_stats["total_bases"]) if read2_stats else 0.0)
    )
    / total_bases
    if total_bases
    else 0.0
)
mean_quality = (
    (
        read1_stats["mean_quality"] * read1_stats["total_bases"]
        + ((read2_stats["mean_quality"] * read2_stats["total_bases"]) if read2_stats else 0.0)
    )
    / total_bases
    if total_bases
    else 0.0
)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "layout",
            "read1_reads",
            "read2_reads",
            "total_reads",
            "read1_bases",
            "read2_bases",
            "total_bases",
            "mean_read_length",
            "min_read_length",
            "max_read_length",
            "gc_percent",
            "n_fraction",
            "mean_quality",
            "pair_count_mismatch",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "sample_id": sample_id,
            "layout": layout,
            "read1_reads": read1_stats["read_count"],
            "read2_reads": read2_stats["read_count"] if read2_stats else 0,
            "total_reads": total_reads,
            "read1_bases": read1_stats["total_bases"],
            "read2_bases": read2_stats["total_bases"] if read2_stats else 0,
            "total_bases": total_bases,
            "mean_read_length": f"{mean_read_length:.2f}",
            "min_read_length": min(
                read1_stats["min_read_length"],
                read2_stats["min_read_length"] if read2_stats else read1_stats["min_read_length"],
            ),
            "max_read_length": max(
                read1_stats["max_read_length"],
                read2_stats["max_read_length"] if read2_stats else read1_stats["max_read_length"],
            ),
            "gc_percent": f"{gc_percent:.4f}",
            "n_fraction": f"{n_fraction:.6f}",
            "mean_quality": f"{mean_quality:.4f}",
            "pair_count_mismatch": int(pair_count_mismatch),
        }
    )
