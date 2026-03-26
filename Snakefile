from pathlib import Path

configfile: "config/config.yaml"


def _strip_known_suffixes(path):
    name = Path(path).name
    for suffix in (".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


GENOME_DIR = Path(config["input_dir"])
GENOME_PATTERNS = ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz")

GENOMES = {}
for pattern in GENOME_PATTERNS:
    for fasta in sorted(GENOME_DIR.glob(pattern)):
        sample = _strip_known_suffixes(fasta)
        GENOMES[sample] = str(fasta)

SAMPLES = sorted(GENOMES)

if not SAMPLES:
    print(f"[kmer-phylo-workflow] No FASTA files found in {GENOME_DIR.resolve()}")


include: "workflow/rules/dataset.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/repeats.smk"


rule all:
    input:
        "results/dataset/genome_dataset.tsv",
        "results/qc/qc_summary.tsv",
        "results/repeats/repeat_annotation_summary.tsv",
        expand("results/qc/{sample}.tsv", sample=SAMPLES),
        expand("results/preprocessed/{sample}.fa", sample=SAMPLES),
        expand("results/repeats/annotation/{sample}.intervals.txt", sample=SAMPLES),
        expand("results/repeats/masked/{sample}.fa", sample=SAMPLES)
