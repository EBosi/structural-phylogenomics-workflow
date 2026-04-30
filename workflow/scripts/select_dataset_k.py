import csv
import math
from pathlib import Path


sample_scores_path = Path(snakemake.input.sample_scores)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

preferred_k_values = sorted({int(value) for value in snakemake.params.preferred_k_values})
min_supported_k = int(snakemake.params.min_supported_k)
support_fraction = float(snakemake.params.support_fraction)
if not 0 < support_fraction <= 1:
    raise ValueError("reads.dataset_k_support_fraction must be in the interval (0, 1]")

rows = list(csv.DictReader(sample_scores_path.open(), delimiter="\t"))
included_rows = [row for row in rows if row["recommended_action"] != "exclude"]
included_count = len(included_rows)

recommended_k = min_supported_k
supporting_samples = 0
required_support = max(1, math.ceil(included_count * support_fraction)) if included_count else 0
selection_status = "no_included_samples" if not included_count else "no_supported_preferred_k"

candidate_k_values = sorted(set(preferred_k_values + [min_supported_k]), reverse=True)
if included_rows:
    for k_value in candidate_k_values:
        support = sum(1 for row in included_rows if int(row["max_supported_k"]) >= k_value)
        if support >= required_support:
            recommended_k = k_value
            supporting_samples = support
            selection_status = "selected"
            break

if not supporting_samples and included_rows:
    supporting_samples = sum(1 for row in included_rows if int(row["max_supported_k"]) >= recommended_k)

dataset_confidence = (
    "high"
    if included_count and supporting_samples == included_count
    else "medium"
    if included_count and supporting_samples >= required_support
    else "low"
)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "n_samples",
            "n_included_samples",
            "n_excluded_samples",
            "recommended_k",
            "required_support",
            "supporting_samples",
            "dataset_confidence",
            "selection_status",
            "selection_method",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "n_samples": len(rows),
            "n_included_samples": included_count,
            "n_excluded_samples": len(rows) - included_count,
            "recommended_k": recommended_k,
            "required_support": required_support,
            "supporting_samples": supporting_samples,
            "dataset_confidence": dataset_confidence,
            "selection_status": selection_status,
            "selection_method": "quality_filtered_histogram_support",
        }
    )
