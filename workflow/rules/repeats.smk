rule repeat_annotation:
    input:
        outpath("organelle", "filtered", "{accession}.fa")
    output:
        outpath("repeats", "annotation", "{accession}.intervals.txt")
    params:
        backend=config["repeat_annotation"]["backend"],
        dustmasker=config["repeat_annotation"]["dustmasker_path"],
        window=config["repeat_annotation"]["dust_window"],
        level=config["repeat_annotation"]["dust_level"],
        linker=config["repeat_annotation"]["dust_linker"],
        repeatmasker=config["repeat_annotation"]["repeatmasker_path"],
        repeatmasker_species=config["repeat_annotation"]["repeatmasker_species"],
        repeatmasker_library=config["repeat_annotation"]["repeatmasker_library"],
        repeatmasker_engine=config["repeat_annotation"]["repeatmasker_engine"],
        repeatmasker_threads=config["repeat_annotation"]["repeatmasker_threads"],
        repeatmasker_extra_args=config["repeat_annotation"]["repeatmasker_extra_args"]
    threads: 1
    script:
        "../scripts/run_repeat_annotation.py"


rule repeat_annotation_sample_summary:
    input:
        fasta=outpath("organelle", "filtered", "{accession}.fa"),
        intervals=outpath("repeats", "annotation", "{accession}.intervals.txt")
    output:
        outpath("repeats", "annotation", "{accession}.summary.tsv")
    params:
        sample="{accession}"
    script:
        "../scripts/summarize_repeat_annotation.py"


rule repeat_annotation_summary:
    input:
        expand(outpath("repeats", "annotation", "{accession}.summary.tsv"), accession=SAMPLE_IDS)
    output:
        outpath("repeats", "repeat_annotation_summary.tsv")
    script:
        "../scripts/merge_tables.py"


rule repeat_masking:
    input:
        fasta=outpath("organelle", "filtered", "{accession}.fa"),
        intervals=outpath("repeats", "annotation", "{accession}.intervals.txt")
    output:
        outpath("repeats", "masked", "{accession}.fa")
    params:
        hard_masking=config["repeat_annotation"]["hard_masking"]
    script:
        "../scripts/mask_fasta_from_intervals.py"
