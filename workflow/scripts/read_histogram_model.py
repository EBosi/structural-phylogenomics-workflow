from __future__ import annotations

from collections import defaultdict
from math import ceil
from typing import Iterable


def group_histogram_rows(rows: Iterable[dict[str, str]]) -> dict[tuple[str, int], list[dict[str, str]]]:
    grouped: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["sample_id"], int(row["k"]))].append(row)
    return dict(grouped)


def summarize_hist_rows(hist_rows: list[dict[str, str]]) -> dict[str, str]:
    if not hist_rows:
        raise ValueError("Cannot summarize an empty histogram row group")

    metadata = hist_rows[0]
    distinct_kmers = int(metadata["distinct_kmers"])
    bins = sorted(
        (int(row["count_bin"]), int(row["n_distinct_kmers_at_count"]))
        for row in hist_rows
        if int(row["count_bin"]) > 0
    )

    if distinct_kmers == 0 or not bins:
        return {
            "sample_id": metadata["sample_id"],
            "k": metadata["k"],
            "sampled_reads": metadata["sampled_reads"],
            "total_kmers": metadata["total_kmers"],
            "distinct_kmers": metadata["distinct_kmers"],
            "singleton_fraction": "0.000000",
            "retained_fraction": "0.000000",
            "error_peak_abundance": 1,
            "error_peak_count_of_counts": 0,
            "signal_peak_abundance": 0,
            "signal_peak_count_of_counts": 0,
            "low_band_suggestion": 2,
            "high_band_suggestion": 3,
            "high_copy_fraction": "0.000000",
            "model_warning": "no_distinct_kmers",
        }

    count_by_abundance = dict(bins)
    singleton_count = count_by_abundance.get(1, 0)
    signal_peak_abundance, signal_peak_count = _find_signal_peak(bins)
    low_band = _find_low_band(bins, signal_peak_abundance)
    high_band = _find_high_band(bins, signal_peak_abundance, signal_peak_count, low_band)

    retained = sum(count for abundance, count in bins if low_band <= abundance <= high_band)
    high_copy = sum(count for abundance, count in bins if abundance > high_band)

    warnings = []
    if signal_peak_abundance == 0:
        warnings.append("no_signal_peak")
    if singleton_count / distinct_kmers > 0.95:
        warnings.append("singleton_dominated")
    if retained / distinct_kmers < 0.05:
        warnings.append("low_retained_fraction")
    if high_copy / distinct_kmers > 0.80:
        warnings.append("high_copy_dominated")

    return {
        "sample_id": metadata["sample_id"],
        "k": metadata["k"],
        "sampled_reads": metadata["sampled_reads"],
        "total_kmers": metadata["total_kmers"],
        "distinct_kmers": metadata["distinct_kmers"],
        "singleton_fraction": f"{singleton_count / distinct_kmers:.6f}",
        "retained_fraction": f"{retained / distinct_kmers:.6f}",
        "error_peak_abundance": 1,
        "error_peak_count_of_counts": singleton_count,
        "signal_peak_abundance": signal_peak_abundance,
        "signal_peak_count_of_counts": signal_peak_count,
        "low_band_suggestion": low_band,
        "high_band_suggestion": high_band,
        "high_copy_fraction": f"{high_copy / distinct_kmers:.6f}",
        "model_warning": ";".join(warnings) if warnings else "ok",
    }


def band_confidence(summary: dict[str, str], recommended_action: str, sample_supports_dataset_k: bool) -> str:
    if recommended_action == "exclude":
        return "excluded"
    if not sample_supports_dataset_k:
        return "unsupported_dataset_k"

    signal_peak = int(summary["signal_peak_abundance"])
    singleton_fraction = float(summary["singleton_fraction"])
    retained_fraction = float(summary["retained_fraction"])
    high_copy_fraction = float(summary["high_copy_fraction"])

    if signal_peak >= 8 and singleton_fraction < 0.70 and retained_fraction >= 0.20 and high_copy_fraction < 0.50:
        return "high"
    if signal_peak >= 3 and singleton_fraction < 0.90 and retained_fraction >= 0.08 and high_copy_fraction < 0.80:
        return "medium"
    return "low"


def _find_signal_peak(bins: list[tuple[int, int]]) -> tuple[int, int]:
    signal_bins = [(abundance, count) for abundance, count in bins if abundance >= 2]
    if not signal_bins:
        return 0, 0

    local_peaks = []
    count_by_abundance = dict(bins)
    for abundance, count in signal_bins:
        previous_count = count_by_abundance.get(abundance - 1, -1)
        next_count = count_by_abundance.get(abundance + 1, -1)
        if count >= previous_count and count >= next_count:
            local_peaks.append((abundance, count))

    candidates = local_peaks or signal_bins
    return max(candidates, key=lambda item: (item[1], -item[0]))


def _find_low_band(bins: list[tuple[int, int]], signal_peak_abundance: int) -> int:
    if signal_peak_abundance <= 2:
        return 2

    valley_candidates = [
        (abundance, count)
        for abundance, count in bins
        if 2 <= abundance <= signal_peak_abundance
    ]
    if not valley_candidates:
        return 2

    valley_abundance, _ = min(valley_candidates, key=lambda item: (item[1], item[0]))
    return max(2, valley_abundance)


def _find_high_band(
    bins: list[tuple[int, int]],
    signal_peak_abundance: int,
    signal_peak_count: int,
    low_band: int,
) -> int:
    if signal_peak_abundance <= 0:
        return max(low_band + 1, 3)

    tail_floor = max(1, ceil(signal_peak_count * 0.01))
    for abundance, count in bins:
        if abundance <= signal_peak_abundance:
            continue
        if count <= tail_floor:
            return max(low_band + 1, abundance)

    highest_observed = max(abundance for abundance, _ in bins)
    empirical_cap = max(low_band + 1, signal_peak_abundance * 3)
    return min(highest_observed, empirical_cap) if highest_observed > low_band else low_band + 1
