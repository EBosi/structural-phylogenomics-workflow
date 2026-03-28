from pathlib import Path

from io_utils import read_fasta


def normalize_header(header):
    return header.split()[0]

def wrap_sequence(sequence, width=80):
    for idx in range(0, len(sequence), width):
        yield sequence[idx : idx + width]


sample = snakemake.params.sample
input_fasta = snakemake.input[0]
output_fasta = Path(snakemake.output[0])
output_fasta.parent.mkdir(parents=True, exist_ok=True)

min_contig_length = int(snakemake.params.min_contig_length)
normalize_headers = bool(snakemake.params.normalize_headers)
uppercase_sequences = bool(snakemake.params.uppercase_sequences)

written = 0

with output_fasta.open("w") as out_handle:
    for index, (header, sequence) in enumerate(read_fasta(input_fasta), start=1):
        if len(sequence) < min_contig_length:
            continue

        clean_header = normalize_header(header) if normalize_headers else header
        clean_sequence = sequence.upper() if uppercase_sequences else sequence

        out_handle.write(f">{sample}|contig_{index}|{clean_header}\n")
        for chunk in wrap_sequence(clean_sequence):
            out_handle.write(f"{chunk}\n")
        written += 1

if written == 0:
    raise ValueError(
        f"No sequences left after preprocessing for sample '{sample}'. "
        "Relax the filters in config/config.yaml."
    )
