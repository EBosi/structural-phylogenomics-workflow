from pathlib import Path

from io_utils import read_fasta


def write_fasta_prefix(input_fasta, output_fasta, max_bases):
    written = 0
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w") as handle:
        for header, sequence in read_fasta(input_fasta):
            if written >= max_bases:
                break
            remaining = max_bases - written
            trimmed = sequence[:remaining]
            if not trimmed:
                continue
            handle.write(f">{header}\n")
            for start in range(0, len(trimmed), 80):
                handle.write(trimmed[start : start + 80] + "\n")
            written += len(trimmed)
    if written == 0:
        raise ValueError(f"No sequence bases available for FASTA prefix subset: {input_fasta}")


write_fasta_prefix(
    snakemake.input.fasta,
    snakemake.output[0],
    int(snakemake.params.max_bases),
)
