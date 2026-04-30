import csv
from pathlib import Path

from read_histogram_model import band_confidence, group_histogram_rows, summarize_hist_rows


histogram_path = Path(snakemake.input.histogram)
model_path = Path(snakemake.input.model)
dataset_k_path = Path(snakemake.input.dataset_k)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

hist_rows = list(csv.DictReader(histogram_path.open(), delimiter="\t"))
model_row = next(csv.DictReader(model_path.open(), delimiter="\t"))
dataset_row = next(csv.DictReader(dataset_k_path.open(), delimiter="\t"))

selected_k = int(dataset_row["recommended_k"])
grouped_hist_rows = group_histogram_rows(hist_rows)
k_rows = grouped_hist_rows.get((model_row["sample_id"], selected_k), [])
if not k_rows:
    raise ValueError(f"No histogram rows found for dataset k={selected_k} in {histogram_path}")

summary = summarize_hist_rows(k_rows)
max_supported_k = int(model_row["max_supported_k"])
recommended_action = model_row["recommended_action"]
sample_supports_dataset_k = recommended_action != "exclude" and max_supported_k >= selected_k
confidence = band_confidence(summary, recommended_action, sample_supports_dataset_k)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "selected_k",
            "low_count",
            "high_count",
            "signal_peak_abundance",
            "signal_peak_count_of_counts",
            "singleton_fraction",
            "retained_fraction",
            "high_copy_fraction",
            "sample_supports_dataset_k",
            "recommended_action",
            "model_warning",
            "band_confidence",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "sample_id": model_row["sample_id"],
            "selected_k": selected_k,
            "low_count": summary["low_band_suggestion"],
            "high_count": summary["high_band_suggestion"],
            "signal_peak_abundance": summary["signal_peak_abundance"],
            "signal_peak_count_of_counts": summary["signal_peak_count_of_counts"],
            "singleton_fraction": summary["singleton_fraction"],
            "retained_fraction": summary["retained_fraction"],
            "high_copy_fraction": summary["high_copy_fraction"],
            "sample_supports_dataset_k": int(sample_supports_dataset_k),
            "recommended_action": recommended_action,
            "model_warning": summary["model_warning"],
            "band_confidence": confidence,
        }
    )
