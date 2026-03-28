import csv
from pathlib import Path


def load_table(path, key):
    rows = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row[key]] = row
    return rows


assemblies = load_table(snakemake.input.assemblies, "accession")
qc = load_table(snakemake.input.qc, "sample")
preprocessing = load_table(snakemake.input.preprocessing, "sample")
organelle = load_table(snakemake.input.organelle, "sample")
repeats = load_table(snakemake.input.repeats, "sample")

rows = []
for accession in sorted(assemblies):
    assembly = assemblies[accession]
    qc_row = qc.get(accession, {})
    prep_row = preprocessing.get(accession, {})
    organelle_row = organelle.get(accession, {})
    repeat_row = repeats.get(accession, {})
    raw_bases = int(qc_row.get("total_bases", 0) or 0)
    raw_sequences = int(qc_row.get("n_sequences", 0) or 0)
    post_preprocess_bases = int(prep_row.get("processed_bases", 0) or 0)
    post_preprocess_sequences = int(prep_row.get("processed_sequences", 0) or 0)
    organelle_removed_bases = int(organelle_row.get("removed_bases", 0) or 0)
    organelle_removed_sequences = int(organelle_row.get("removed_sequences", 0) or 0)
    post_organelle_bases = int(organelle_row.get("kept_bases", 0) or 0)
    post_organelle_sequences = int(organelle_row.get("kept_sequences", 0) or 0)
    post_organelle_base_fraction = round((post_organelle_bases / raw_bases), 6) if raw_bases else 0
    post_organelle_sequence_fraction = round((post_organelle_sequences / raw_sequences), 6) if raw_sequences else 0
    rows.append(
        {
            "accession": accession,
            "resolved_accession": assembly.get("resolved_accession", ""),
            "organism_name": assembly.get("organism_name", ""),
            "source_db": assembly.get("source_db", ""),
            "assembly_level": assembly.get("assembly_level", ""),
            "raw_sequences": raw_sequences,
            "raw_bases": raw_bases,
            "raw_n50": qc_row.get("n50", ""),
            "raw_gc_percent": qc_row.get("gc_percent", ""),
            "post_preprocess_sequences": post_preprocess_sequences,
            "post_preprocess_bases": post_preprocess_bases,
            "post_preprocess_sequence_fraction": prep_row.get("retained_sequence_fraction", ""),
            "post_preprocess_base_fraction": prep_row.get("retained_base_fraction", ""),
            "post_organelle_sequences": post_organelle_sequences,
            "post_organelle_bases": post_organelle_bases,
            "post_organelle_sequence_fraction": post_organelle_sequence_fraction,
            "post_organelle_base_fraction": post_organelle_base_fraction,
            "organelle_removed_sequences": organelle_removed_sequences,
            "organelle_removed_bases": organelle_removed_bases,
            "organelle_confident_sequences": organelle_row.get("organelle_confident_sequences", ""),
            "organelle_ambiguous_sequences": organelle_row.get("organelle_ambiguous_sequences", ""),
            "masked_intervals": repeat_row.get("masked_intervals", ""),
            "masked_bases": repeat_row.get("masked_bases", ""),
            "masked_fraction_percent": repeat_row.get("masked_fraction_percent", ""),
        }
    )

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
