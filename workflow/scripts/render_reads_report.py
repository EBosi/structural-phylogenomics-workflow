import csv
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_signature_status(path: Path) -> dict[str, str]:
    rows = read_table(path)
    if not rows:
        return {"sample_id": path.stem, "signature_status": "missing_or_empty", "hashes": "0"}
    sample_id = rows[0].get("accession", path.stem)
    status = rows[0].get("signature_status", "ok")
    hashes = sum(1 for row in rows if row.get("hash_value"))
    retained = rows[0].get("retained_kmers", "0")
    return {
        "sample_id": sample_id,
        "signature_status": status,
        "hashes": str(hashes),
        "retained_kmers": retained,
    }


def distance_matrix_size(path: Path) -> int:
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    return max(0, len(header) - 1)


manifest_rows = read_table(Path(snakemake.input.manifest))
qc_rows = read_table(Path(snakemake.input.qc_summary))
model_rows = read_table(Path(snakemake.input.model_summary))
histogram_rows = read_table(Path(snakemake.input.histogram_summary))
score_rows = read_table(Path(snakemake.input.sample_scores))
dataset_row = read_table(Path(snakemake.input.dataset_k))[0]
band_rows = [next(csv.DictReader(Path(path).open(), delimiter="\t")) for path in snakemake.input.band_reports]
signature_rows = [read_signature_status(Path(path)) for path in snakemake.input.signature_reports]
sketch_distance_path = Path(snakemake.input.sketch_distance)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

sketch_sample_count = distance_matrix_size(sketch_distance_path)
recommended_k_note = dataset_row["recommended_k"]
if dataset_row.get("selection_status") == "no_included_samples":
    recommended_k_note = f"{recommended_k_note} (fallback; no included samples)"

lines = [
    "# Reads-Mode Report",
    "",
    "This target is an experimental adaptive read-based phylogenomics workflow.",
    "It performs manifest validation, FASTQ QC, sampled k-mer count histograms, histogram-guided sample scoring, dataset-level K selection, abundance-band filtering, filtered read sketching, and safe tree inference for retained samples.",
    "",
    "## Dataset Summary",
    "",
    f"- Samples in manifest: {len(manifest_rows)}",
    f"- Included samples before sketching: {dataset_row['n_included_samples']}",
    f"- Excluded samples before sketching: {dataset_row['n_excluded_samples']}",
    f"- Samples retained in sketch distance matrix: {sketch_sample_count}",
    f"- Recommended dataset K: {recommended_k_note}",
    f"- Dataset confidence: {dataset_row['dataset_confidence']}",
    f"- K selection status: {dataset_row.get('selection_status', '')}",
    "",
    "## Current Limitations",
    "",
    "- Histogram modeling is peak-aware but still heuristic; it is not a fitted mixture/coverage model.",
    "- Histograms and filtered sketches are computed from bounded read prefixes, so file ordering can bias estimates.",
    "- Python `Counter`-based k-mer counting is not suitable for large publication-scale read sets without an external counter or streaming sketch/count-min design.",
    "- Samples with failed QC, unsupported dataset K, low-confidence bands, or zero retained k-mers are removed from the sketch distance matrix.",
    "",
    "## Sample Overview",
    "",
]

score_by_sample = {row["sample_id"]: row for row in score_rows}
band_by_sample = {row["sample_id"]: row for row in band_rows}
signature_by_sample = {row["sample_id"]: row for row in signature_rows}
hist_by_sample = {}
for row in histogram_rows:
    sample_id = row["sample_id"]
    current = hist_by_sample.get(sample_id)
    if current is None or int(row["k"]) > int(current["k"]):
        hist_by_sample[sample_id] = row
for qc_row in qc_rows:
    sample_id = qc_row["sample_id"]
    score_row = score_by_sample.get(sample_id, {})
    hist_row = hist_by_sample.get(sample_id, {})
    band_row = band_by_sample.get(sample_id, {})
    signature_row = signature_by_sample.get(sample_id, {})
    lines.append(
        f"- {sample_id}: reads={qc_row['total_reads']}, bases={qc_row['total_bases']}, "
        f"mean_len={qc_row['mean_read_length']}, mean_q={qc_row.get('mean_quality', 'NA')}, "
        f"N_frac={qc_row.get('n_fraction', 'NA')}, coverage={score_row.get('approx_coverage', '0')}, "
        f"max_k={score_row.get('max_supported_k', '')}, signal_peak={hist_row.get('signal_peak_abundance', '0')}, "
        f"band={band_row.get('low_count', '0')}-{band_row.get('high_count', '0')}, "
        f"band_conf={band_row.get('band_confidence', '')}, action={score_row.get('recommended_action', '')}, "
        f"signature={signature_row.get('signature_status', 'missing')}, hashes={signature_row.get('hashes', '0')}"
    )

lines.append("")
lines.append("## Notes")
lines.append("")
lines.append("- Coverage is estimated from total read bases and user-provided estimated genome size.")
lines.append("- Dataset K is selected from samples that pass QC; per-sample bands are inferred at that dataset K.")
lines.append("- This is still not publication-grade until validated on real FASTQ sets with explicit sensitivity analyses against coverage, read length, repeats, contamination, and known phylogenies.")

output_path.write_text("\n".join(lines) + "\n")
