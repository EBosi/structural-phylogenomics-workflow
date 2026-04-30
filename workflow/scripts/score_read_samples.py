import csv
from pathlib import Path


manifest_path = Path(snakemake.input.manifest)
model_summary_path = Path(snakemake.input.model_summary)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

manifest_rows = {}
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        manifest_rows[row["sample_id"]] = row

rows = []
with model_summary_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        sample_id = row["sample_id"]
        manifest_row = manifest_rows.get(sample_id, {})
        rows.append(
            {
                "sample_id": sample_id,
                "species": manifest_row.get("species", ""),
                "approx_coverage": row["approx_coverage"],
                "quality_band": row["quality_band"],
                "max_supported_k": row["max_supported_k"],
                "signal_peak_abundance": row.get("signal_peak_abundance", ""),
                "retained_fraction": row.get("retained_fraction", ""),
                "singleton_fraction": row.get("singleton_fraction", ""),
                "high_copy_fraction": row.get("high_copy_fraction", ""),
                "sample_quality_score": row["sample_quality_score"],
                "recommended_action": row["recommended_action"],
                "failure_reasons": row.get("failure_reasons", ""),
                "platform": manifest_row.get("platform", ""),
            }
        )

rows.sort(key=lambda row: row["sample_id"])

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "species",
            "platform",
            "approx_coverage",
            "quality_band",
            "max_supported_k",
            "signal_peak_abundance",
            "retained_fraction",
            "singleton_fraction",
            "high_copy_fraction",
            "sample_quality_score",
            "recommended_action",
            "failure_reasons",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
