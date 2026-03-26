rule qc_per_sample:
    input:
        "data/genomes/{accession}.fna.gz"
    output:
        "results/qc/{accession}.tsv"
    params:
        sample="{accession}"
    script:
        "../scripts/qc_stats.py"


rule qc_summary:
    input:
        expand("results/qc/{accession}.tsv", accession=ACCESSIONS)
    output:
        "results/qc/qc_summary.tsv"
    script:
        "../scripts/merge_tables.py"
