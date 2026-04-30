import csv
from pathlib import Path


def jaccard_distance(set_a, set_b):
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return 1.0 - (len(set_a & set_b) / union)


def input_signature_paths():
    if hasattr(snakemake.input, "signatures"):
        return list(snakemake.input.signatures)
    return list(snakemake.input)


def requested_names():
    names = getattr(snakemake.params, "accessions", None)
    if names is None:
        names = getattr(snakemake.params, "samples", [])
    return list(names)


signatures = {}
signature_paths = input_signature_paths()
requested = requested_names()
excluded_statuses = set(getattr(snakemake.params, "exclude_signature_statuses", []))

for index, path in enumerate(signature_paths):
    expected_name = requested[index] if index < len(requested) else None
    current_accession = None
    values = set()
    statuses = set()
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            current_accession = row.get("accession") or current_accession
            if row.get("signature_status"):
                statuses.add(row["signature_status"])
            hash_value = row.get("hash_value", "")
            if hash_value != "":
                values.add(int(hash_value))

    accession = current_accession or expected_name
    if accession is None:
        continue
    if statuses & excluded_statuses:
        continue
    signatures[accession] = values

accessions = [accession for accession in requested if accession in signatures]
if not accessions and not requested:
    accessions = sorted(signatures)

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession"] + accessions)
    for acc_a in accessions:
        row = [acc_a]
        for acc_b in accessions:
            row.append(f"{jaccard_distance(signatures[acc_a], signatures[acc_b]):.12f}")
        writer.writerow(row)
