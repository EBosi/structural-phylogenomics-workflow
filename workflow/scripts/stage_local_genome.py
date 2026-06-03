import gzip
import shutil
from pathlib import Path


input_fasta = Path(snakemake.input[0])
output_fasta = Path(snakemake.output[0])
output_fasta.parent.mkdir(parents=True, exist_ok=True)

if input_fasta.resolve() == output_fasta.resolve():
    raise ValueError(f"Input and output FASTA paths are identical: {input_fasta}")

if str(input_fasta).lower().endswith(".gz"):
    with input_fasta.open("rb") as in_handle, output_fasta.open("wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)
else:
    with input_fasta.open("rb") as in_handle, gzip.open(output_fasta, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)
