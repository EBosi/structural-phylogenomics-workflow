import csv
import subprocess
from pathlib import Path

from io_utils import read_fasta


def fasta_lengths(path):
    return {header.split()[0]: len(sequence) for header, sequence in read_fasta(path)}


query_path = Path(snakemake.input.fasta)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

tmp_path = output_path.with_suffix(".blast.tsv")
blast_cmd = [
    snakemake.params.blastn,
    "-query",
    str(query_path),
    "-db",
    snakemake.params.db_prefix,
    "-outfmt",
    "6 qseqid sseqid pident length evalue bitscore qlen",
    "-max_target_seqs",
    "5",
    "-evalue",
    "1e-10",
    "-task",
    "blastn",
    "-out",
    str(tmp_path),
]
subprocess.run(blast_cmd, check=True)

lengths = fasta_lengths(query_path)
best_hits = {}

if tmp_path.exists():
    with tmp_path.open() as handle:
        for raw_line in handle:
            qseqid, sseqid, pident, length, evalue, bitscore, qlen = raw_line.rstrip("\n").split("\t")
            pident = float(pident)
            length = int(length)
            evalue = float(evalue)
            bitscore = float(bitscore)
            qlen = int(qlen)
            qcov = (length / qlen) * 100 if qlen else 0.0
            current = best_hits.get(qseqid)
            candidate = {
                "query_id": qseqid,
                "query_length": qlen,
                "subject_id": sseqid,
                "pident": pident,
                "aligned_length": length,
                "query_coverage_percent": qcov,
                "evalue": evalue,
                "bitscore": bitscore,
            }
            if current is None or (candidate["bitscore"], candidate["aligned_length"]) > (
                current["bitscore"],
                current["aligned_length"],
            ):
                best_hits[qseqid] = candidate

rows = []
for query_id, qlen in sorted(lengths.items()):
    hit = best_hits.get(query_id)
    if hit is None:
        rows.append(
            {
                "query_id": query_id,
                "query_length": qlen,
                "subject_id": "",
                "pident": "",
                "aligned_length": 0,
                "query_coverage_percent": 0,
                "evalue": "",
                "bitscore": "",
                "classification": "nuclear_like",
            }
        )
        continue

    if (
        hit["pident"] >= float(snakemake.params.confident_identity)
        and hit["query_coverage_percent"] >= float(snakemake.params.confident_qcov)
    ):
        classification = "organelle_confident"
    elif (
        hit["pident"] >= float(snakemake.params.ambiguous_identity)
        and hit["query_coverage_percent"] >= float(snakemake.params.ambiguous_qcov)
    ):
        classification = "organelle_ambiguous"
    else:
        classification = "nuclear_like"

    rows.append(
        {
            "query_id": query_id,
            "query_length": qlen,
            "subject_id": hit["subject_id"],
            "pident": f"{hit['pident']:.4f}",
            "aligned_length": hit["aligned_length"],
            "query_coverage_percent": f"{hit['query_coverage_percent']:.4f}",
            "evalue": hit["evalue"],
            "bitscore": f"{hit['bitscore']:.4f}",
            "classification": classification,
        }
    )

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "query_id",
            "query_length",
            "subject_id",
            "pident",
            "aligned_length",
            "query_coverage_percent",
            "evalue",
            "bitscore",
            "classification",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
