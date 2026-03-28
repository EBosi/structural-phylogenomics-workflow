import csv
from collections import Counter
from pathlib import Path

from io_utils import read_fasta


def wrap_sequence(sequence, width=80):
    for idx in range(0, len(sequence), width):
        yield sequence[idx : idx + width]


remove_classes = set(snakemake.params.remove_classes)
accession = snakemake.params.accession

classifications = {}
with Path(snakemake.input.calls).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        classifications[row["query_id"]] = row["classification"]

kept_sequences = 0
kept_bases = 0
removed_sequences = 0
removed_bases = 0
class_counter = Counter(classifications.values())

output_fasta = Path(snakemake.output.filtered)
output_fasta.parent.mkdir(parents=True, exist_ok=True)

with output_fasta.open("w") as out_handle:
    for header, sequence in read_fasta(snakemake.input.fasta):
        query_id = header.split()[0]
        classification = classifications.get(query_id, "nuclear_like")
        if classification in remove_classes:
            removed_sequences += 1
            removed_bases += len(sequence)
            continue
        kept_sequences += 1
        kept_bases += len(sequence)
        out_handle.write(f">{header}\n")
        for chunk in wrap_sequence(sequence):
            out_handle.write(f"{chunk}\n")

summary = {
    "sample": accession,
    "kept_sequences": kept_sequences,
    "removed_sequences": removed_sequences,
    "kept_bases": kept_bases,
    "removed_bases": removed_bases,
    "nuclear_like_sequences": class_counter.get("nuclear_like", 0),
    "organelle_ambiguous_sequences": class_counter.get("organelle_ambiguous", 0),
    "organelle_confident_sequences": class_counter.get("organelle_confident", 0),
}

summary_path = Path(snakemake.output.summary)
summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary.keys()), delimiter="\t")
    writer.writeheader()
    writer.writerow(summary)
