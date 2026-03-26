rule dataset_table:
    output:
        "results/dataset/genome_dataset.tsv"
    params:
        genomes=GENOMES
    script:
        "../scripts/build_dataset_table.py"
