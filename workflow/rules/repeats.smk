rule repeat_annotation:
    input:
        "results/preprocessed/{sample}.fa"
    output:
        "results/repeats/annotation/{sample}.intervals.txt"
    params:
        dustmasker=config["repeat_annotation"]["dustmasker_path"],
        window=config["repeat_annotation"]["dust_window"],
        level=config["repeat_annotation"]["dust_level"],
        linker=config["repeat_annotation"]["dust_linker"]
    shell:
        (
            "{params.dustmasker} "
            "-in {input} "
            "-out {output} "
            "-outfmt interval "
            "-window {params.window} "
            "-level {params.level} "
            "-linker {params.linker}"
        )


rule repeat_annotation_sample_summary:
    input:
        fasta="results/preprocessed/{sample}.fa",
        intervals="results/repeats/annotation/{sample}.intervals.txt"
    output:
        "results/repeats/annotation/{sample}.summary.tsv"
    params:
        sample="{sample}"
    script:
        "../scripts/summarize_repeat_annotation.py"


rule repeat_annotation_summary:
    input:
        expand("results/repeats/annotation/{sample}.summary.tsv", sample=SAMPLES)
    output:
        "results/repeats/repeat_annotation_summary.tsv"
    script:
        "../scripts/merge_tables.py"


rule repeat_masking:
    input:
        "results/preprocessed/{sample}.fa"
    output:
        "results/repeats/masked/{sample}.fa"
    params:
        dustmasker=config["repeat_annotation"]["dustmasker_path"],
        window=config["repeat_annotation"]["dust_window"],
        level=config["repeat_annotation"]["dust_level"],
        linker=config["repeat_annotation"]["dust_linker"],
        hard_masking=config["repeat_annotation"]["hard_masking"],
        hard_masking_flag="-hard_masking" if config["repeat_annotation"]["hard_masking"] else ""
    shell:
        (
            "{params.dustmasker} "
            "-in {input} "
            "-out {output} "
            "-outfmt fasta "
            "-window {params.window} "
            "-level {params.level} "
            "-linker {params.linker} "
            "{params.hard_masking_flag}"
        )
