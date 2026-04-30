import csv
from collections import defaultdict
from pathlib import Path

from io_utils import read_fasta
from repeat_utils import add_interval, load_interval_details, merge_intervals


SOURCE_COLUMNS = [
    "dustmasker",
    "repeatmasker_known",
    "repeatmasker_custom",
    "repeatmasker_denovo",
]


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


def grouped_interval_stats(records, key_fields):
    grouped = defaultdict(lambda: defaultdict(list))
    for record in records:
        key = tuple(record.get(field, "") for field in key_fields)
        add_interval(grouped[key], record["seq_id"], int(record["start"]), int(record["end"]))

    stats = {}
    for key, intervals_by_seq in grouped.items():
        merged = []
        for intervals in intervals_by_seq.values():
            merged.extend(merge_intervals(intervals))
        stats[key] = {
            "intervals": len(merged),
            "bases": sum(end - start + 1 for start, end in merged),
        }
    return stats


sample = snakemake.params.sample
seq_lengths = load_sequence_lengths(snakemake.input.fasta)
interval_count, masked_bases = parse_interval_file(snakemake.input.intervals)
total_bases = sum(seq_lengths.values())
records = load_interval_details(snakemake.input.details)
source_stats = grouped_interval_stats(records, ["source"])
class_stats = grouped_interval_stats(records, ["source", "repeat_class", "repeat_family"])

row = {
    "sample": sample,
    "n_sequences": len(seq_lengths),
    "total_bases": total_bases,
    "masked_intervals": interval_count,
    "masked_bases": masked_bases,
    "masked_fraction_percent": round((masked_bases / total_bases) * 100, 4) if total_bases else 0,
}

for source in SOURCE_COLUMNS:
    stats = source_stats.get((source,), {"intervals": 0, "bases": 0})
    row[f"{source}_intervals"] = stats["intervals"]
    row[f"{source}_bases"] = stats["bases"]
    row[f"{source}_fraction_percent"] = round((stats["bases"] / total_bases) * 100, 4) if total_bases else 0

class_rows = []
for (source, repeat_class, repeat_family), stats in sorted(class_stats.items()):
    class_rows.append(
        {
            "sample": sample,
            "source": source,
            "repeat_class": repeat_class or "Unknown",
            "repeat_family": repeat_family,
            "masked_intervals": stats["intervals"],
            "masked_bases": stats["bases"],
            "masked_fraction_percent": round((stats["bases"] / total_bases) * 100, 4) if total_bases else 0,
        }
    )

summary_path = Path(snakemake.output.summary)
summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(row)

classes_path = Path(snakemake.output.classes)
classes_path.parent.mkdir(parents=True, exist_ok=True)
with classes_path.open("w", newline="") as handle:
    fieldnames = [
        "sample",
        "source",
        "repeat_class",
        "repeat_family",
        "masked_intervals",
        "masked_bases",
        "masked_fraction_percent",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(class_rows)
