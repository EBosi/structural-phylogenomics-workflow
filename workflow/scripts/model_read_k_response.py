import csv
from pathlib import Path

from read_histogram_model import group_histogram_rows, summarize_hist_rows


def parse_float(value: str, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except ValueError:
        return default


def parse_bool_int(value: str) -> bool:
    return str(value).strip() in {"1", "true", "True", "yes", "YES"}


def coverage_limited_k(
    candidate_ks: list[int],
    approx_coverage: float,
    low_cov: float,
    medium_cov: float,
    high_cov: float,
    warn_cov: float,
) -> int:
    ordered = sorted(candidate_ks)
    if approx_coverage >= high_cov:
        return ordered[-1]
    if approx_coverage >= medium_cov:
        return ordered[max(0, len(ordered) - 2)]
    if approx_coverage >= low_cov:
        return ordered[min(len(ordered) - 1, 1)]
    if approx_coverage >= warn_cov:
        return ordered[0]
    return ordered[0]


def choose_fallback_k(
    candidate_ks: list[int],
    mean_read_length: float,
    approx_coverage: float,
    min_mean_read_length: int,
    low_cov: float,
    medium_cov: float,
    high_cov: float,
    warn_cov: float,
) -> int:
    coverage_cap = coverage_limited_k(candidate_ks, approx_coverage, low_cov, medium_cov, high_cov, warn_cov)
    admissible = [
        k_value
        for k_value in candidate_ks
        if k_value <= coverage_cap and mean_read_length >= max(min_mean_read_length, 2 * k_value)
    ]
    return max(admissible) if admissible else min(candidate_ks)


manifest = Path(snakemake.input.manifest)
qc_path = Path(snakemake.input.qc)
histogram_path = Path(snakemake.input.histogram)
sample_id = snakemake.wildcards.sample
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

candidate_k_values = sorted(int(value) for value in snakemake.params.candidate_k_values)
min_supported_k = int(snakemake.params.min_supported_k)
low_cov = float(snakemake.params.low_coverage_threshold)
medium_cov = float(snakemake.params.medium_coverage_threshold)
high_cov = float(snakemake.params.high_coverage_threshold)
min_mean_read_length = int(snakemake.params.min_mean_read_length)
warn_cov = float(snakemake.params.include_with_caution_coverage_threshold)
min_retained_fraction = float(getattr(snakemake.params, "min_retained_fraction", 0.05))
max_singleton_fraction = float(getattr(snakemake.params, "max_singleton_fraction", 0.98))
max_high_copy_fraction = float(getattr(snakemake.params, "max_high_copy_fraction", 0.85))
min_mean_quality = float(getattr(snakemake.params, "min_mean_quality", 0.0))
max_n_fraction = float(getattr(snakemake.params, "max_n_fraction", 1.0))
exclude_pair_mismatch = bool(getattr(snakemake.params, "exclude_pair_mismatch", True))

if min_supported_k not in candidate_k_values:
    raise ValueError("reads.min_supported_k must be present in reads.candidate_k_values")

manifest_row = None
with manifest.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["sample_id"] == sample_id:
            manifest_row = row
            break

if manifest_row is None:
    raise ValueError(f"Sample {sample_id} not found in manifest {manifest}")

with qc_path.open() as handle:
    qc_row = next(csv.DictReader(handle, delimiter="\t"))

estimated_genome_size_bp = parse_float(manifest_row["estimated_genome_size_bp"])
total_bases = parse_float(qc_row["total_bases"])
mean_read_length = parse_float(qc_row["mean_read_length"])
mean_quality = parse_float(qc_row.get("mean_quality", "0"))
n_fraction = parse_float(qc_row.get("n_fraction", "0"))
pair_count_mismatch = parse_bool_int(qc_row.get("pair_count_mismatch", "0"))
approx_coverage = (total_bases / estimated_genome_size_bp) if estimated_genome_size_bp > 0 else 0.0

all_hist_rows = list(csv.DictReader(histogram_path.open(), delimiter="\t"))
grouped_hist_rows = group_histogram_rows(all_hist_rows)
hist_summaries = []
for k_value in candidate_k_values:
    k_rows = grouped_hist_rows.get((sample_id, k_value), [])
    if k_rows:
        hist_summaries.append(summarize_hist_rows(k_rows))

supported_k_values = []
selected_summary = None
for summary in hist_summaries:
    k_value = int(summary["k"])
    singleton_fraction = float(summary["singleton_fraction"])
    retained_fraction = float(summary["retained_fraction"])
    signal_peak_abundance = int(summary["signal_peak_abundance"])
    high_copy_fraction = float(summary["high_copy_fraction"])
    if mean_read_length < max(min_mean_read_length, 2 * k_value):
        continue
    if retained_fraction < min_retained_fraction:
        continue
    if singleton_fraction > max_singleton_fraction:
        continue
    if high_copy_fraction > max_high_copy_fraction:
        continue
    if approx_coverage >= warn_cov and signal_peak_abundance >= 2:
        supported_k_values.append(k_value)

if supported_k_values:
    max_supported_k = max(supported_k_values)
else:
    max_supported_k = max(
        min_supported_k,
        choose_fallback_k(
            candidate_k_values,
            mean_read_length,
            approx_coverage,
            min_mean_read_length,
            low_cov,
            medium_cov,
            high_cov,
            warn_cov,
        ),
    )
    max_supported_k = max(k for k in candidate_k_values if k <= max_supported_k)

for summary in hist_summaries:
    if int(summary["k"]) == max_supported_k:
        selected_summary = summary
        break

if selected_summary is None and hist_summaries:
    selected_summary = hist_summaries[0]
    max_supported_k = int(selected_summary["k"])

if approx_coverage >= high_cov:
    quality_band = "high"
elif approx_coverage >= medium_cov:
    quality_band = "medium"
elif approx_coverage >= low_cov:
    quality_band = "low"
else:
    quality_band = "very_low"

failure_reasons = []
if estimated_genome_size_bp <= 0:
    failure_reasons.append("invalid_genome_size")
if mean_read_length < min_mean_read_length:
    failure_reasons.append("short_reads")
if approx_coverage < warn_cov:
    failure_reasons.append("coverage_below_include_threshold")
if exclude_pair_mismatch and pair_count_mismatch:
    failure_reasons.append("paired_read_count_mismatch")
if min_mean_quality > 0 and mean_quality < min_mean_quality:
    failure_reasons.append("low_mean_quality")
if n_fraction > max_n_fraction:
    failure_reasons.append("high_n_fraction")
if hist_summaries and not supported_k_values and approx_coverage >= warn_cov:
    failure_reasons.append("no_histogram_supported_k")

if failure_reasons:
    recommended_action = "exclude"
elif approx_coverage < low_cov:
    recommended_action = "include_with_caution"
else:
    recommended_action = "include"

signal_peak_abundance = int(selected_summary["signal_peak_abundance"]) if selected_summary else 0
singleton_fraction = float(selected_summary["singleton_fraction"]) if selected_summary else 1.0
retained_fraction = float(selected_summary["retained_fraction"]) if selected_summary else 0.0
high_copy_fraction = float(selected_summary["high_copy_fraction"]) if selected_summary else 1.0
model_warning = selected_summary["model_warning"] if selected_summary else "no_histogram_summary"

informative_low_count = int(selected_summary["low_band_suggestion"]) if selected_summary else 2
informative_high_count = int(selected_summary["high_band_suggestion"]) if selected_summary else 3
coverage_score = min(1.0, max(0.0, approx_coverage / high_cov if high_cov else 0.0))
histogram_score = max(0.0, min(1.0, retained_fraction * (1.0 - singleton_fraction) * (1.0 - high_copy_fraction)))
read_quality_score = min(1.0, max(0.0, mean_quality / 30.0)) if mean_quality else 0.0
sample_quality_score = (0.55 * coverage_score) + (0.35 * histogram_score) + (0.10 * read_quality_score)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "estimated_genome_size_bp",
            "approx_coverage",
            "quality_band",
            "max_supported_k",
            "informative_low_count",
            "informative_high_count",
            "signal_peak_abundance",
            "singleton_fraction",
            "retained_fraction",
            "high_copy_fraction",
            "mean_quality",
            "n_fraction",
            "pair_count_mismatch",
            "sample_quality_score",
            "recommended_action",
            "failure_reasons",
            "model_warning",
            "model_type",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "sample_id": sample_id,
            "estimated_genome_size_bp": int(estimated_genome_size_bp) if estimated_genome_size_bp else 0,
            "approx_coverage": f"{approx_coverage:.4f}",
            "quality_band": quality_band,
            "max_supported_k": max_supported_k,
            "informative_low_count": informative_low_count,
            "informative_high_count": informative_high_count,
            "signal_peak_abundance": signal_peak_abundance,
            "singleton_fraction": f"{singleton_fraction:.6f}",
            "retained_fraction": f"{retained_fraction:.6f}",
            "high_copy_fraction": f"{high_copy_fraction:.6f}",
            "mean_quality": f"{mean_quality:.4f}",
            "n_fraction": f"{n_fraction:.6f}",
            "pair_count_mismatch": int(pair_count_mismatch),
            "sample_quality_score": f"{sample_quality_score:.4f}",
            "recommended_action": recommended_action,
            "failure_reasons": ";".join(failure_reasons) if failure_reasons else "none",
            "model_warning": model_warning,
            "model_type": "shared_peak_aware_histogram_guided_v2",
        }
    )
