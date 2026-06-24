import csv
import re
from pathlib import Path

configfile: "config/config.yaml"
OUTPUT_ROOT = str(config.get("output_root", "results")).rstrip("/")
GENOME_ROOT = str(config.get("genome_root", "data/genomes")).rstrip("/")
RESOURCE_ROOT = str(config.get("resource_root", "resources")).rstrip("/")


def outpath(*parts):
    return str(Path(OUTPUT_ROOT, *parts))


def genome_path(filename):
    return str(Path(GENOME_ROOT, filename))


def resource_path(*parts):
    return str(Path(RESOURCE_ROOT, *parts))


ACCESSION_FILE = Path(config["metadata"]["accession_file"])
ACCESSION_CLI = config.get("accessions", "")
LOCAL_GENOME_FILE = Path(config["metadata"]["local_genomes_file"])
LOCAL_GENOMES_DIR = config.get("local_genomes_dir") or config["metadata"].get("local_genomes_dir", "")
LOCAL_GENOME_EXTENSIONS = tuple(
    config["metadata"].get(
        "local_genome_extensions",
        [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"],
    )
)


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


def _matches_local_genome_suffix(path):
    name = path.name.lower()
    return any(name.endswith(suffix.lower()) for suffix in LOCAL_GENOME_EXTENSIONS)


def _sample_id_from_local_genome_path(path):
    name = path.name
    for suffix in sorted(LOCAL_GENOME_EXTENSIONS, key=len, reverse=True):
        if name.lower().endswith(suffix.lower()):
            sample_id = name[: -len(suffix)]
            break
    else:
        raise ValueError(f"Local genome file {path} does not match configured FASTA extensions")

    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_id):
        raise ValueError(
            f"Local genome file {path} yields invalid sample ID '{sample_id}'. "
            "Use filenames containing only letters, numbers, dots, underscores and hyphens."
        )
    return sample_id


def _discover_local_genomes_dir():
    if not LOCAL_GENOMES_DIR:
        return {}

    local_dir = Path(LOCAL_GENOMES_DIR)
    if not local_dir.exists():
        raise ValueError(f"Configured local_genomes_dir does not exist: {local_dir}")
    if not local_dir.is_dir():
        raise ValueError(f"Configured local_genomes_dir is not a directory: {local_dir}")

    discovered = {}
    for fasta_path in sorted(local_dir.iterdir()):
        if not fasta_path.is_file() or not _matches_local_genome_suffix(fasta_path):
            continue
        sample_id = _sample_id_from_local_genome_path(fasta_path)
        if sample_id in discovered:
            raise ValueError(
                f"Multiple local genome files in {local_dir} resolve to sample ID '{sample_id}'"
            )
        discovered[sample_id] = str(fasta_path)

    if not discovered:
        suffixes = ", ".join(LOCAL_GENOME_EXTENSIONS)
        raise ValueError(f"No local genome FASTA files found in {local_dir}; expected suffixes: {suffixes}")

    return discovered


def _combine_local_sample_ids(tsv_sample_ids, dir_sample_ids):
    overlap = sorted(set(tsv_sample_ids) & set(dir_sample_ids))
    if overlap:
        overlap_str = ", ".join(overlap)
        raise ValueError(
            "The following sample identifiers are present in both metadata/local_genomes.tsv "
            f"and local_genomes_dir: {overlap_str}"
        )
    return list(tsv_sample_ids) + list(dir_sample_ids)


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
LOCAL_TSV_SAMPLE_IDS = _parse_local_sample_ids()
LOCAL_DIR_FASTAS = _discover_local_genomes_dir()
LOCAL_DIR_SAMPLE_IDS = list(LOCAL_DIR_FASTAS)
LOCAL_SAMPLE_IDS = _combine_local_sample_ids(LOCAL_TSV_SAMPLE_IDS, LOCAL_DIR_SAMPLE_IDS)
SAMPLE_IDS = _combine_sample_ids(REMOTE_ACCESSIONS, LOCAL_SAMPLE_IDS)
LOCAL_DIR_SAMPLE_PATTERN = "|".join(re.escape(accession) for accession in LOCAL_DIR_SAMPLE_IDS) or r"(?!)"
LOCAL_DIR_RECORDS = [
    {
        "accession": accession,
        "organism_name": accession,
        "assembly_name": accession,
        "assembly_level": "local",
        "source_db": "local",
        "local_path": genome_path(f"{accession}.fna.gz"),
    }
    for accession in LOCAL_DIR_SAMPLE_IDS
]
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
        "Populate metadata/accessions.txt, metadata/local_genomes.tsv, "
        "or pass --config accessions=... or local_genomes_dir=..."
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
        outpath("reports", "pre_kmer_report.md"),
        outpath("reports", "pre_kmer_summary.tsv")


rule pre_kmer:
    input:
        outpath("metadata", "assemblies.tsv"),
        outpath("metadata", "organisms.tsv"),
        outpath("metadata", "download_manifest.tsv"),
        outpath("qc", "qc_summary.tsv"),
        outpath("preprocessing", "preprocessing_summary.tsv"),
        outpath("organelle", "organelle_summary.tsv"),
        outpath("repeats", "repeat_annotation_summary.tsv"),
        outpath("reports", "pre_kmer_report.md"),
        outpath("reports", "pre_kmer_summary.tsv"),
        expand(genome_path("{accession}.fna.gz"), accession=SAMPLE_IDS),
        expand(outpath("qc", "{accession}.tsv"), accession=SAMPLE_IDS),
        expand(outpath("preprocessed", "{accession}.fa"), accession=SAMPLE_IDS),
        expand(outpath("organelle", "calls", "{accession}.tsv"), accession=SAMPLE_IDS),
        expand(outpath("organelle", "filtered", "{accession}.fa"), accession=SAMPLE_IDS),
        expand(outpath("organelle", "{accession}.summary.tsv"), accession=SAMPLE_IDS),
        expand(outpath("repeats", "annotation", "{accession}.intervals.txt"), accession=SAMPLE_IDS),
        expand(outpath("repeats", "annotation", "{accession}.summary.tsv"), accession=SAMPLE_IDS),
        expand(outpath("preprocessing", "{accession}.summary.tsv"), accession=SAMPLE_IDS),
        expand(outpath("repeats", "masked", "{accession}.fa"), accession=SAMPLE_IDS)


rule full_analysis:
    input:
        rules.pre_kmer.input,
        expand(
            outpath("kmers", "matrices", "{dataset}", "k{k}.tsv"),
            dataset=KMER_DATASETS,
            k=K_VALUES,
        ),
        expand(
            outpath("distances", "{dataset}", "k{k}", "{metric}.tsv"),
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
        ),
        expand(
            outpath("trees", "{dataset}", "k{k}", "{metric}", "{method}.nwk"),
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
            method=TREE_METHODS,
        ),
        outpath("reports", "tree_manifest.tsv"),
        outpath("reports", "tree_comparisons.tsv"),
        expand(
            outpath("resampling", "bootstrap", "{dataset}", "k{k}", "{metric}", "{method}.summary.tsv"),
            dataset=RESAMPLING_DATASETS,
            k=RESAMPLING_K_VALUES,
            metric=RESAMPLING_METRICS,
            method=RESAMPLING_METHODS,
        ),
        expand(
            outpath("resampling", "jackknife", "{dataset}", "k{k}", "{metric}", "{method}.summary.tsv"),
            dataset=RESAMPLING_DATASETS,
            k=RESAMPLING_K_VALUES,
            metric=RESAMPLING_METRICS,
            method=RESAMPLING_METHODS,
        ),
        expand(
            outpath("sketch", "distances", "{dataset}", "k{k}", "minhash_jaccard.tsv"),
            dataset=SKETCH_DATASETS,
            k=SKETCH_K_VALUES,
        ),
        expand(
            outpath("sketch", "trees", "{dataset}", "k{k}", "{method}.nwk"),
            dataset=SKETCH_DATASETS,
            k=SKETCH_K_VALUES,
            method=SKETCH_METHODS,
        )
