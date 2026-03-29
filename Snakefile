import csv
from pathlib import Path

configfile: "config/config.yaml"
ACCESSION_FILE = Path(config["metadata"]["accession_file"])
ACCESSION_CLI = config.get("accessions", "")
LOCAL_GENOME_FILE = Path(config["metadata"]["local_genomes_file"])


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


def _parse_local_sample_ids():
    if not LOCAL_GENOME_FILE.exists():
        return []

    sample_ids = []
    with LOCAL_GENOME_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        if "accession" not in reader.fieldnames:
            raise ValueError(
                f"Local genome metadata file {LOCAL_GENOME_FILE} must contain an 'accession' column"
            )
        for row in reader:
            accession = row["accession"].strip()
            if accession:
                sample_ids.append(accession)

    deduplicated = []
    seen = set()
    for accession in sample_ids:
        if accession not in seen:
            deduplicated.append(accession)
            seen.add(accession)
    return deduplicated


def _combine_sample_ids(remote_accessions, local_sample_ids):
    overlap = sorted(set(remote_accessions) & set(local_sample_ids))
    if overlap:
        overlap_str = ", ".join(overlap)
        raise ValueError(
            "The following sample identifiers are present in both remote accessions and local genomes: "
            f"{overlap_str}"
        )

    combined = []
    seen = set()
    for accession in list(remote_accessions) + list(local_sample_ids):
        if accession not in seen:
            combined.append(accession)
            seen.add(accession)
    return combined


REMOTE_ACCESSIONS = _parse_accessions()
LOCAL_SAMPLE_IDS = _parse_local_sample_ids()
SAMPLE_IDS = _combine_sample_ids(REMOTE_ACCESSIONS, LOCAL_SAMPLE_IDS)
K_VALUES = [int(value) for value in config["kmers"]["k_values"]]
KMER_DATASETS = list(config["kmers"]["datasets"])
DISTANCE_METRICS = list(config["distances"]["metrics"])
TREE_METHODS = list(config["trees"]["methods"])
RESAMPLING_DATASETS = list(config["resampling"]["datasets"])
RESAMPLING_K_VALUES = [int(value) for value in config["resampling"]["k_values"]]
RESAMPLING_METRICS = list(config["resampling"]["metrics"])
RESAMPLING_METHODS = list(config["resampling"]["methods"])
SKETCH_DATASETS = list(config["sketch"]["datasets"])
SKETCH_K_VALUES = [int(value) for value in config["sketch"]["k_values"]]
SKETCH_METHODS = list(config["sketch"]["methods"])

if not SAMPLE_IDS:
    print(
        "[kmer-phylo-workflow] No samples found. "
        "Populate metadata/accessions.txt, metadata/local_genomes.tsv, or pass --config accessions=..."
    )


include: "workflow/rules/metadata.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/organelle.smk"
include: "workflow/rules/repeats.smk"
include: "workflow/rules/pre_kmer_reports.smk"
include: "workflow/rules/kmers.smk"
include: "workflow/rules/distances.smk"
include: "workflow/rules/trees.smk"
include: "workflow/rules/reports.smk"
include: "workflow/rules/resampling.smk"
include: "workflow/rules/sketch.smk"


rule all:
    input:
        "results/reports/pre_kmer_report.md",
        "results/reports/pre_kmer_summary.tsv"


rule pre_kmer:
    input:
        "results/metadata/assemblies.tsv",
        "results/metadata/organisms.tsv",
        "results/metadata/download_manifest.tsv",
        "results/qc/qc_summary.tsv",
        "results/preprocessing/preprocessing_summary.tsv",
        "results/organelle/organelle_summary.tsv",
        "results/repeats/repeat_annotation_summary.tsv",
        "results/reports/pre_kmer_report.md",
        "results/reports/pre_kmer_summary.tsv",
        expand("data/genomes/{accession}.fna.gz", accession=SAMPLE_IDS),
        expand("results/qc/{accession}.tsv", accession=SAMPLE_IDS),
        expand("results/preprocessed/{accession}.fa", accession=SAMPLE_IDS),
        expand("results/organelle/calls/{accession}.tsv", accession=SAMPLE_IDS),
        expand("results/organelle/filtered/{accession}.fa", accession=SAMPLE_IDS),
        expand("results/organelle/{accession}.summary.tsv", accession=SAMPLE_IDS),
        expand("results/repeats/annotation/{accession}.intervals.txt", accession=SAMPLE_IDS),
        expand("results/repeats/annotation/{accession}.summary.tsv", accession=SAMPLE_IDS),
        expand("results/preprocessing/{accession}.summary.tsv", accession=SAMPLE_IDS),
        expand("results/repeats/masked/{accession}.fa", accession=SAMPLE_IDS)


rule full_analysis:
    input:
        rules.pre_kmer.input,
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
        "results/reports/tree_comparisons.tsv",
        expand(
            "results/resampling/bootstrap/{dataset}/k{k}/{metric}/{method}.summary.tsv",
            dataset=RESAMPLING_DATASETS,
            k=RESAMPLING_K_VALUES,
            metric=RESAMPLING_METRICS,
            method=RESAMPLING_METHODS,
        ),
        expand(
            "results/resampling/jackknife/{dataset}/k{k}/{metric}/{method}.summary.tsv",
            dataset=RESAMPLING_DATASETS,
            k=RESAMPLING_K_VALUES,
            metric=RESAMPLING_METRICS,
            method=RESAMPLING_METHODS,
        ),
        expand(
            "results/sketch/distances/{dataset}/k{k}/minhash_jaccard.tsv",
            dataset=SKETCH_DATASETS,
            k=SKETCH_K_VALUES,
        ),
        expand(
            "results/sketch/trees/{dataset}/k{k}/{method}.nwk",
            dataset=SKETCH_DATASETS,
            k=SKETCH_K_VALUES,
            method=SKETCH_METHODS,
        )
