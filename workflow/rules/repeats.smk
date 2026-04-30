def repeat_backend_components():
    return {
        component.strip().lower()
        for component in str(config["repeat_annotation"]["backend"]).split("+")
        if component.strip()
    }


def repeatmasker_pilot_accessions():
    values = config["repeat_annotation"].get("repeatmasker_pilot_accessions", [])
    if isinstance(values, str):
        return [value.strip() for value in values.split(",") if value.strip()]
    return [str(value).strip() for value in values if str(value).strip()]


def load_repeatmasker_sample_config():
    import csv
    from pathlib import Path

    path = Path(config["repeat_annotation"].get("repeatmasker_sample_config", ""))
    if not path or not path.exists():
        return {}

    rows = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        if "accession" not in reader.fieldnames:
            raise ValueError(f"RepeatMasker sample config {path} must contain an accession column")
        for row in reader:
            accession = row.get("accession", "").strip()
            if accession:
                rows[accession] = {key: (value.strip() if value is not None else "") for key, value in row.items()}
    return rows


REPEATMASKER_SAMPLE_CONFIG = load_repeatmasker_sample_config()
REPEATMASKER_PILOT_ACCESSIONS = repeatmasker_pilot_accessions()


def repeatmasker_setting(wildcards, key):
    sample_row = REPEATMASKER_SAMPLE_CONFIG.get(wildcards.accession, {})
    value = sample_row.get(key, "")
    if value:
        return value
    return config["repeat_annotation"].get(key, "")


def repeatmodeler_library_input(wildcards):
    if "repeatmodeler" in repeat_backend_components():
        return f"results/repeats/repeatmodeler/{wildcards.accession}/consensi.fa.classified"
    return []


rule repeatmodeler_library:
    input:
        fasta="results/organelle/filtered/{accession}.fa"
    output:
        library="results/repeats/repeatmodeler/{accession}/consensi.fa.classified",
        summary="results/repeats/repeatmodeler/{accession}/summary.tsv",
        workdir=directory("results/repeats/repeatmodeler/{accession}/work")
    params:
        repeatmodeler=config["repeat_annotation"]["repeatmodeler_path"],
        build_database=config["repeat_annotation"]["build_database_path"],
        database_name="{accession}",
        ltr_struct=config["repeat_annotation"]["repeatmodeler_ltr_struct"],
        extra_args=config["repeat_annotation"]["repeatmodeler_extra_args"]
    threads: config["repeat_annotation"]["repeatmodeler_threads"]
    script:
        "../scripts/run_repeatmodeler.py"


rule repeat_annotation:
    input:
        fasta="results/organelle/filtered/{accession}.fa",
        denovo_library=repeatmodeler_library_input
    output:
        intervals="results/repeats/annotation/{accession}.intervals.txt",
        details="results/repeats/annotation/{accession}.details.tsv",
        raw_dir=directory("results/repeats/raw/{accession}")
    params:
        backend=config["repeat_annotation"]["backend"],
        dustmasker=config["repeat_annotation"]["dustmasker_path"],
        window=config["repeat_annotation"]["dust_window"],
        level=config["repeat_annotation"]["dust_level"],
        linker=config["repeat_annotation"]["dust_linker"],
        repeatmasker=config["repeat_annotation"]["repeatmasker_path"],
        repeatmasker_species=lambda wc: repeatmasker_setting(wc, "repeatmasker_species"),
        repeatmasker_library=lambda wc: repeatmasker_setting(wc, "repeatmasker_library"),
        repeatmasker_engine=config["repeat_annotation"]["repeatmasker_engine"],
        repeatmasker_threads=config["repeat_annotation"]["repeatmasker_threads"],
        repeatmasker_extra_args=lambda wc: repeatmasker_setting(wc, "repeatmasker_extra_args"),
        repeatmasker_require_species_or_library=config["repeat_annotation"]["repeatmasker_require_species_or_library"],
        max_input_bases=0
    threads: config["repeat_annotation"]["repeatmasker_threads"]
    script:
        "../scripts/run_repeat_annotation.py"


rule repeatmasker_pilot_input:
    input:
        fasta="results/organelle/filtered/{accession}.fa"
    output:
        "results/repeats/pilot/{accession}/input.first{max_bases}.fa"
    params:
        max_bases="{max_bases}"
    script:
        "../scripts/subset_fasta_prefix.py"


rule repeatmasker_pilot_annotation:
    input:
        fasta="results/repeats/pilot/{accession}/input.first{max_bases}.fa"
    output:
        intervals="results/repeats/pilot/{accession}/first{max_bases}/annotation.intervals.txt",
        details="results/repeats/pilot/{accession}/first{max_bases}/annotation.details.tsv",
        raw_dir=directory("results/repeats/pilot/{accession}/first{max_bases}/raw")
    params:
        backend=config["repeat_annotation"].get("repeatmasker_pilot_backend", "dustmasker+repeatmasker"),
        dustmasker=config["repeat_annotation"]["dustmasker_path"],
        window=config["repeat_annotation"]["dust_window"],
        level=config["repeat_annotation"]["dust_level"],
        linker=config["repeat_annotation"]["dust_linker"],
        repeatmasker=config["repeat_annotation"]["repeatmasker_path"],
        repeatmasker_species=lambda wc: repeatmasker_setting(wc, "repeatmasker_species"),
        repeatmasker_library=lambda wc: repeatmasker_setting(wc, "repeatmasker_library"),
        repeatmasker_engine=config["repeat_annotation"]["repeatmasker_engine"],
        repeatmasker_threads=config["repeat_annotation"]["repeatmasker_threads"],
        repeatmasker_extra_args=lambda wc: repeatmasker_setting(wc, "repeatmasker_extra_args"),
        repeatmasker_require_species_or_library=config["repeat_annotation"]["repeatmasker_require_species_or_library"],
        max_input_bases=0
    threads: config["repeat_annotation"]["repeatmasker_threads"]
    script:
        "../scripts/run_repeat_annotation.py"


rule repeatmasker_pilot_sample_summary:
    input:
        fasta="results/repeats/pilot/{accession}/input.first{max_bases}.fa",
        intervals="results/repeats/pilot/{accession}/first{max_bases}/annotation.intervals.txt",
        details="results/repeats/pilot/{accession}/first{max_bases}/annotation.details.tsv"
    output:
        summary="results/repeats/pilot/{accession}/first{max_bases}/summary.tsv",
        classes="results/repeats/pilot/{accession}/first{max_bases}/classes.tsv"
    params:
        sample="{accession}"
    script:
        "../scripts/summarize_repeat_annotation.py"


rule repeatmasker_pilot_summary:
    input:
        expand(
            "results/repeats/pilot/{accession}/first{max_bases}/summary.tsv",
            accession=REPEATMASKER_PILOT_ACCESSIONS,
            max_bases=[config["repeat_annotation"].get("repeatmasker_pilot_max_bases", 1000000)],
        )
    output:
        "results/repeats/pilot/summary.tsv"
    script:
        "../scripts/merge_tables.py"


rule repeatmasker_pilot_class_summary:
    input:
        expand(
            "results/repeats/pilot/{accession}/first{max_bases}/classes.tsv",
            accession=REPEATMASKER_PILOT_ACCESSIONS,
            max_bases=[config["repeat_annotation"].get("repeatmasker_pilot_max_bases", 1000000)],
        )
    output:
        "results/repeats/pilot/class_summary.tsv"
    script:
        "../scripts/merge_tables.py"


rule repeat_annotation_sample_summary:
    input:
        fasta="results/organelle/filtered/{accession}.fa",
        intervals="results/repeats/annotation/{accession}.intervals.txt",
        details="results/repeats/annotation/{accession}.details.tsv"
    output:
        summary="results/repeats/annotation/{accession}.summary.tsv",
        classes="results/repeats/annotation/{accession}.classes.tsv"
    params:
        sample="{accession}"
    script:
        "../scripts/summarize_repeat_annotation.py"


rule repeat_annotation_summary:
    input:
        expand("results/repeats/annotation/{accession}.summary.tsv", accession=SAMPLE_IDS)
    output:
        "results/repeats/repeat_annotation_summary.tsv"
    script:
        "../scripts/merge_tables.py"


rule repeat_class_summary:
    input:
        expand("results/repeats/annotation/{accession}.classes.tsv", accession=SAMPLE_IDS)
    output:
        "results/repeats/repeat_class_summary.tsv"
    script:
        "../scripts/merge_tables.py"


rule repeat_masking:
    input:
        fasta="results/organelle/filtered/{accession}.fa",
        intervals="results/repeats/annotation/{accession}.intervals.txt"
    output:
        "results/repeats/masked/{accession}.fa"
    params:
        hard_masking=config["repeat_annotation"]["hard_masking"]
    script:
        "../scripts/mask_fasta_from_intervals.py"
