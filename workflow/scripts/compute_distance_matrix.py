import csv
import math
from pathlib import Path


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


metric = snakemake.params.metric
input_path = Path(snakemake.input[0])
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with input_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    feature_names = [name for name in reader.fieldnames if name != "accession"]
    data = [(row["accession"], [float(row[name]) for name in feature_names]) for row in reader]

if metric == "cosine":
    distance_fn = cosine_distance
elif metric == "jensen_shannon":
    distance_fn = jensen_shannon_distance
else:
    raise ValueError(f"Unsupported distance metric: {metric}")

accessions = [item[0] for item in data]
vectors = [item[1] for item in data]

with output_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession"] + accessions)
    for accession_a, vector_a in zip(accessions, vectors):
        row = [accession_a]
        for vector_b in vectors:
            row.append(f"{distance_fn(vector_a, vector_b):.12f}")
        writer.writerow(row)
