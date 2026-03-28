rule fetch_organelle_references:
    input:
        "results/metadata/assemblies.tsv"
    output:
        fasta="resources/organelle/mitochondrion_refs.fasta",
        table="resources/organelle/mitochondrion_refs.tsv"
    params:
        esearch=config["organelle_screen"]["esearch_path"],
        efetch=config["organelle_screen"]["efetch_path"],
        organelle_types=config["organelle_screen"]["organelle_types"]
    script:
        "../scripts/fetch_organelle_references.py"


rule build_organelle_blastdb:
    input:
        fasta="resources/organelle/mitochondrion_refs.fasta"
    output:
        nhr="resources/organelle/blastdb/mitochondrion.nhr",
        nin="resources/organelle/blastdb/mitochondrion.nin",
        nsq="resources/organelle/blastdb/mitochondrion.nsq"
    params:
        makeblastdb=config["organelle_screen"]["makeblastdb_path"]
    shell:
        "{params.makeblastdb} -in {input.fasta} -dbtype nucl -out resources/organelle/blastdb/mitochondrion"


rule organelle_screen:
    input:
        fasta="results/preprocessed/{accession}.fa",
        db="resources/organelle/blastdb/mitochondrion.nsq"
    output:
        "results/organelle/calls/{accession}.tsv"
    params:
        blastn=config["organelle_screen"]["blastn_path"],
        db_prefix="resources/organelle/blastdb/mitochondrion",
        confident_identity=config["organelle_screen"]["confident_identity"],
        confident_qcov=config["organelle_screen"]["confident_query_coverage"],
        ambiguous_identity=config["organelle_screen"]["ambiguous_identity"],
        ambiguous_qcov=config["organelle_screen"]["ambiguous_query_coverage"]
    script:
        "../scripts/screen_organelles.py"


rule organelle_filter:
    input:
        fasta="results/preprocessed/{accession}.fa",
        calls="results/organelle/calls/{accession}.tsv"
    output:
        filtered="results/organelle/filtered/{accession}.fa",
        summary="results/organelle/{accession}.summary.tsv"
    params:
        accession="{accession}",
        remove_classes=config["organelle_screen"]["remove_classes"]
    script:
        "../scripts/filter_organelles.py"


rule organelle_summary:
    input:
        expand("results/organelle/{accession}.summary.tsv", accession=ACCESSIONS)
    output:
        "results/organelle/organelle_summary.tsv"
    script:
        "../scripts/merge_tables.py"
