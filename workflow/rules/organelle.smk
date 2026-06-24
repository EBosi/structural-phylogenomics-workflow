rule fetch_organelle_references:
    input:
        outpath("metadata", "assemblies.tsv")
    output:
        fasta=resource_path("organelle", "mitochondrion_refs.fasta"),
        table=resource_path("organelle", "mitochondrion_refs.tsv")
    params:
        esearch=config["organelle_screen"]["esearch_path"],
        efetch=config["organelle_screen"]["efetch_path"],
        organelle_types=config["organelle_screen"]["organelle_types"]
    script:
        "../scripts/fetch_organelle_references.py"


rule build_organelle_blastdb:
    input:
        fasta=resource_path("organelle", "mitochondrion_refs.fasta")
    output:
        nhr=resource_path("organelle", "blastdb", "mitochondrion.nhr"),
        nin=resource_path("organelle", "blastdb", "mitochondrion.nin"),
        nsq=resource_path("organelle", "blastdb", "mitochondrion.nsq")
    params:
        makeblastdb=config["organelle_screen"]["makeblastdb_path"],
        db_prefix=resource_path("organelle", "blastdb", "mitochondrion")
    shell:
        "{params.makeblastdb} -in {input.fasta} -dbtype nucl -out {params.db_prefix}"


rule organelle_screen:
    input:
        fasta=outpath("preprocessed", "{accession}.fa"),
        db=resource_path("organelle", "blastdb", "mitochondrion.nsq")
    output:
        outpath("organelle", "calls", "{accession}.tsv")
    params:
        blastn=config["organelle_screen"]["blastn_path"],
        db_prefix=resource_path("organelle", "blastdb", "mitochondrion"),
        confident_identity=config["organelle_screen"]["confident_identity"],
        confident_qcov=config["organelle_screen"]["confident_query_coverage"],
        ambiguous_identity=config["organelle_screen"]["ambiguous_identity"],
        ambiguous_qcov=config["organelle_screen"]["ambiguous_query_coverage"]
    script:
        "../scripts/screen_organelles.py"


rule organelle_filter:
    input:
        fasta=outpath("preprocessed", "{accession}.fa"),
        calls=outpath("organelle", "calls", "{accession}.tsv")
    output:
        filtered=outpath("organelle", "filtered", "{accession}.fa"),
        summary=outpath("organelle", "{accession}.summary.tsv")
    params:
        accession="{accession}",
        remove_classes=config["organelle_screen"]["remove_classes"]
    script:
        "../scripts/filter_organelles.py"


rule organelle_summary:
    input:
        expand(outpath("organelle", "{accession}.summary.tsv"), accession=SAMPLE_IDS)
    output:
        outpath("organelle", "organelle_summary.tsv")
    script:
        "../scripts/merge_tables.py"
