import csv
from pathlib import Path

from read_histogram_model import group_histogram_rows, summarize_hist_rows


output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

all_rows = []
for histogram_report in snakemake.input.histogram_reports:
    with Path(histogram_report).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        all_rows.extend(reader)

grouped_rows = group_histogram_rows(all_rows)
rows = [summarize_hist_rows(hist_rows) for _, hist_rows in sorted(grouped_rows.items(), key=lambda item: (item[0][0], int(item[0][1])))]

with output_path.open("w", newline="") as handle:
    if rows:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    else:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "k"])
