import csv
from pathlib import Path


assemblies = []
with Path(snakemake.input.assemblies).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    assemblies.extend(reader)

summary_rows = []
with Path(snakemake.input.summary).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    summary_rows.extend(reader)

repeat_class_rows = []
with Path(snakemake.input.repeat_classes).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    repeat_class_rows.extend(reader)

repeat_backend = getattr(snakemake.params, "repeat_backend", "dustmasker")

total_samples = len(summary_rows)
total_raw_bases = sum(int(row["raw_bases"]) for row in summary_rows if row["raw_bases"])
total_post_preprocess_bases = sum(
    int(row["post_preprocess_bases"]) for row in summary_rows if row["post_preprocess_bases"]
)
total_post_organelle_bases = sum(
    int(row["post_organelle_bases"]) for row in summary_rows if row["post_organelle_bases"]
)
total_organelle_removed_bases = sum(
    int(row["organelle_removed_bases"]) for row in summary_rows if row["organelle_removed_bases"]
)

lines = []
lines.append("# Pre-kmer Report")
lines.append("")
lines.append("## Overview")
lines.append("")
lines.append(f"- Samples: {total_samples}")
lines.append(f"- Total raw bases: {total_raw_bases}")
lines.append(f"- Total bases after preprocessing: {total_post_preprocess_bases}")
lines.append(f"- Total bases after organelle filtering: {total_post_organelle_bases}")
lines.append(f"- Total bases removed as organellar: {total_organelle_removed_bases}")
lines.append("")
lines.append("## Sample Summary")
lines.append("")
lines.append(
    "| Accession | Resolved | Organism | Source | Assembly | Raw bases | Post-preprocess bases | Post-organelle bases | Organelle removed bases | Post-organelle retained fraction | Combined masked % | DUST % | RM known % | RM custom % | RM de novo % |"
)
lines.append("| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
for row in summary_rows:
    lines.append(
        "| {accession} | {resolved_accession} | {organism_name} | {source_db} | {assembly_level} | {raw_bases} | {post_preprocess_bases} | {post_organelle_bases} | {organelle_removed_bases} | {post_organelle_base_fraction} | {masked_fraction_percent} | {dustmasker_fraction_percent} | {repeatmasker_known_fraction_percent} | {repeatmasker_custom_fraction_percent} | {repeatmasker_denovo_fraction_percent} |".format(
            **row
        )
    )
lines.append("")
lines.append("## Repeat Class Breakdown")
lines.append("")
if repeat_class_rows:
    lines.append("| Sample | Source | Class | Family | Masked bases | Masked % |")
    lines.append("| --- | --- | --- | --- | ---: | ---: |")
    for row in repeat_class_rows:
        lines.append(
            "| {sample} | {source} | {repeat_class} | {repeat_family} | {masked_bases} | {masked_fraction_percent} |".format(
                **row
            )
        )
else:
    lines.append("No repeat class annotations were produced for this run.")
lines.append("")
lines.append("## Notes")
lines.append("")
lines.append("- Raw QC is computed on input genome FASTA files, whether downloaded or supplied locally.")
lines.append("- Preprocessing removes short contigs and normalizes FASTA formatting.")
lines.append("- Organelle screening uses BLAST against mitochondrial reference sequences fetched from NCBI.")
lines.append("- Organelle filtering removes only contigs classified as `organelle_confident`.")
lines.append(f"- Repeat annotation/masking backend for this run: `{repeat_backend}`.")
lines.append("- `dustmasker` is only a low-complexity/simple-sequence filter; it is not TE annotation.")
lines.append("- RepeatMasker class labels are retained when available in `results/repeats/repeat_class_summary.tsv`.")
lines.append("- Combined masked percentage is the union of all interval sources, so it can be lower than the sum of source-specific percentages when sources overlap.")

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text("\n".join(lines) + "\n")
