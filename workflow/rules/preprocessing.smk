rule preprocess_genome:
    input:
        "data/genomes/{accession}.fna.gz"
    output:
        "results/preprocessed/{accession}.fa"
    params:
        sample="{accession}",
        min_contig_length=config["preprocessing"]["min_contig_length"],
        normalize_headers=config["preprocessing"]["normalize_headers"],
        uppercase_sequences=config["preprocessing"]["uppercase_sequences"],
        exclude_keywords=config["preprocessing"]["exclude_organelle_keywords"]
    script:
        "../scripts/preprocess_genome.py"
