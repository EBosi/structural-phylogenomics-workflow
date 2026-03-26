rule qc_per_sample:
    input:
        lambda wc: GENOMES[wc.sample]
    output:
        "results/qc/{sample}.tsv"
    params:
        sample="{sample}"
    script:
        "../scripts/qc_stats.py"


rule qc_summary:
    input:
        expand("results/qc/{sample}.tsv", sample=SAMPLES)
    output:
        "results/qc/qc_summary.tsv"
    script:
        "../scripts/merge_tables.py"
