import csv
import itertools
from pathlib import Path

from Bio import Phylo


def terminal_names(tree):
    return sorted(clade.name for clade in tree.get_terminals())


def split_set(tree):
    taxa = set(terminal_names(tree))
    splits = set()
    for clade in tree.get_nonterminals(order="postorder"):
        leaves = {leaf.name for leaf in clade.get_terminals()}
        if not leaves or leaves == taxa:
            continue
        complement = taxa - leaves
        if len(leaves) < 2 or len(complement) < 2:
            continue
        smaller = tuple(sorted(leaves if len(leaves) <= len(complement) else complement))
        splits.add(smaller)
    return splits


manifest_rows = []
with Path(snakemake.input.manifest).open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    manifest_rows.extend(reader)

loaded = []
for row in manifest_rows:
    tree = Phylo.read(row["tree_path"], "newick")
    loaded.append(
        {
            **row,
            "taxa": terminal_names(tree),
            "splits": split_set(tree),
        }
    )

comparison_rows = []
for left, right in itertools.combinations(loaded, 2):
    same_taxa = left["taxa"] == right["taxa"]
    if same_taxa:
        shared = len(left["splits"] & right["splits"])
        union = len(left["splits"] | right["splits"])
        jaccard = (shared / union) if union else 1.0
    else:
        shared = ""
        union = ""
        jaccard = ""

    comparison_rows.append(
        {
            "left_dataset": left["dataset"],
            "left_k": left["k"],
            "left_metric": left["metric"],
            "left_method": left["method"],
            "right_dataset": right["dataset"],
            "right_k": right["k"],
            "right_metric": right["metric"],
            "right_method": right["method"],
            "same_taxa": same_taxa,
            "shared_splits": shared,
            "union_splits": union,
            "split_jaccard": f"{jaccard:.12f}" if jaccard != "" else "",
        }
    )

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "left_dataset",
            "left_k",
            "left_metric",
            "left_method",
            "right_dataset",
            "right_k",
            "right_metric",
            "right_method",
            "same_taxa",
            "shared_splits",
            "union_splits",
            "split_jaccard",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(comparison_rows)
