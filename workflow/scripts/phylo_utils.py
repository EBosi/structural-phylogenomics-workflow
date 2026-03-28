import csv
import hashlib
import math
from collections import Counter

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def canonical_kmer(kmer):
    rc = reverse_complement(kmer)
    return min(kmer, rc)


def cosine_distance(vec_a, vec_b):
    dot = sum(a * b for a, b in zip(vec_a, vec_b))
    norm_a = math.sqrt(sum(a * a for a in vec_a))
    norm_b = math.sqrt(sum(b * b for b in vec_b))
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return 1.0 - (dot / (norm_a * norm_b))


def jensen_shannon_distance(vec_a, vec_b):
    sum_a = sum(vec_a)
    sum_b = sum(vec_b)
    if sum_a == 0 or sum_b == 0:
        return 0.0

    p = [value / sum_a for value in vec_a]
    q = [value / sum_b for value in vec_b]
    m = [(pa + qa) / 2.0 for pa, qa in zip(p, q)]

    def kl_divergence(x, y):
        total = 0.0
        for xi, yi in zip(x, y):
            if xi == 0:
                continue
            total += xi * math.log2(xi / yi)
        return total

    js_div = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
    return math.sqrt(max(js_div, 0.0))


def distance_function(metric):
    if metric == "cosine":
        return cosine_distance
    if metric == "jensen_shannon":
        return jensen_shannon_distance
    raise ValueError(f"Unsupported distance metric: {metric}")


def build_lower_triangle(names, full_matrix):
    lower = []
    for i, _ in enumerate(names):
        lower.append(full_matrix[i][: i + 1])
    return lower


def infer_tree_from_distance_matrix(names, full_matrix, method):
    distance_matrix = DistanceMatrix(names=names, matrix=build_lower_triangle(names, full_matrix))
    constructor = DistanceTreeConstructor()
    if method == "nj":
        return constructor.nj(distance_matrix)
    if method == "upgma":
        return constructor.upgma(distance_matrix)
    raise ValueError(f"Unsupported tree method: {method}")


def tree_terminal_names(tree):
    return sorted(clade.name for clade in tree.get_terminals())


def tree_split_set(tree):
    taxa = set(tree_terminal_names(tree))
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


def compute_full_distance_matrix(vectors, metric):
    dist_fn = distance_function(metric)
    matrix = []
    for vec_a in vectors:
        row = []
        for vec_b in vectors:
            row.append(dist_fn(vec_a, vec_b))
        matrix.append(row)
    return matrix


def sha1_u64(text):
    digest = hashlib.sha1(text.encode("ascii")).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def aggregate_kmer_rows(rows):
    counts = Counter()
    for row in rows:
        counts[row["kmer"]] += float(row["frequency"])
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: value / total for kmer, value in counts.items()}
