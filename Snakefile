from pathlib import Path

configfile: "config/config.yaml"
ACCESSION_FILE = Path(config["metadata"]["accession_file"])
ACCESSION_CLI = config.get("accessions", "")


def _parse_accessions():
    accession_values = []

    if ACCESSION_FILE.exists():
        accession_values.extend(
            line.strip() for line in ACCESSION_FILE.read_text().splitlines() if line.strip()
        )

    if ACCESSION_CLI:
        accession_values.extend(
            token.strip() for token in ACCESSION_CLI.split(",") if token.strip()
        )

    deduplicated = []
    seen = set()
    for accession in accession_values:
        if accession not in seen:
            deduplicated.append(accession)
            seen.add(accession)
    return deduplicated


ACCESSIONS = _parse_accessions()
K_VALUES = [int(value) for value in config["kmers"]["k_values"]]
KMER_DATASETS = list(config["kmers"]["datasets"])
DISTANCE_METRICS = list(config["distances"]["metrics"])
TREE_METHODS = list(config["trees"]["methods"])

if not ACCESSIONS:
    print(
        "[kmer-phylo-workflow] No accessions found. "
        "Populate metadata/accessions.txt or pass --config accessions=..."
    )


include: "workflow/rules/metadata.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/repeats.smk"
include: "workflow/rules/kmers.smk"
include: "workflow/rules/distances.smk"
include: "workflow/rules/trees.smk"
include: "workflow/rules/reports.smk"


rule all:
    input:
        "results/metadata/assemblies.tsv",
        "results/metadata/organisms.tsv",
        "results/metadata/download_manifest.tsv",
        "results/qc/qc_summary.tsv",
        "results/repeats/repeat_annotation_summary.tsv",
        expand("data/genomes/{accession}.fna.gz", accession=ACCESSIONS),
        expand("results/qc/{accession}.tsv", accession=ACCESSIONS),
        expand("results/preprocessed/{accession}.fa", accession=ACCESSIONS),
        expand("results/repeats/annotation/{accession}.intervals.txt", accession=ACCESSIONS),
        expand("results/repeats/masked/{accession}.fa", accession=ACCESSIONS),
        expand(
            "results/kmers/matrices/{dataset}/k{k}.tsv",
            dataset=KMER_DATASETS,
            k=K_VALUES,
        ),
        expand(
            "results/distances/{dataset}/k{k}/{metric}.tsv",
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
        ),
        expand(
            "results/trees/{dataset}/k{k}/{metric}/{method}.nwk",
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
            method=TREE_METHODS,
        ),
        "results/reports/tree_manifest.tsv",
        "results/reports/tree_comparisons.tsv"
