from repeat_utils import load_interval_text, write_masked_fasta


intervals_by_seq = load_interval_text(snakemake.input.intervals)
write_masked_fasta(
    snakemake.input.fasta,
    snakemake.output[0],
    intervals_by_seq,
    bool(snakemake.params.hard_masking),
)
