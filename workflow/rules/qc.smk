rule qc_per_sample:
    wildcard_constraints:
        accession="(?!qc_summary\\.tsv$)[^/]+"
    input:
        genome_path("{accession}.fna.gz")
    output:
        outpath("qc", "{accession}.tsv")
    params:
        sample="{accession}"
    script:
        "../scripts/qc_stats.py"


rule qc_summary:
    input:
        expand(outpath("qc", "{accession}.tsv"), accession=SAMPLE_IDS)
    output:
        outpath("qc", "qc_summary.tsv")
    script:
        "../scripts/merge_tables.py"
