import csv
from pathlib import Path


def parse_tree_path(path):
    parts = Path(path).parts
    try:
        trees_idx = parts.index("trees")
    except ValueError as exc:
        raise ValueError(f"Could not parse tree path: {path}") from exc

    dataset = parts[trees_idx + 1]
    k_value = parts[trees_idx + 2]
    metric = parts[trees_idx + 3]
    filename = parts[trees_idx + 4]
    method = filename.replace(".nwk", "")
    return dataset, k_value, metric, method


rows = []
for tree_path in sorted(snakemake.input):
    dataset, k_value, metric, method = parse_tree_path(tree_path)
    rows.append(
        {
            "dataset": dataset,
            "k": k_value.replace("k", ""),
            "metric": metric,
            "method": method,
            "tree_path": tree_path,
        }
    )

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["dataset", "k", "metric", "method", "tree_path"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
