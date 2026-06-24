rule preprocess_genome:
    input:
        genome_path("{accession}.fna.gz")
    output:
        outpath("preprocessed", "{accession}.fa")
    params:
        sample="{accession}",
        min_contig_length=config["preprocessing"]["min_contig_length"],
        normalize_headers=config["preprocessing"]["normalize_headers"],
        uppercase_sequences=config["preprocessing"]["uppercase_sequences"]
    script:
        "../scripts/preprocess_genome.py"


rule preprocess_summary_per_sample:
    input:
        raw=genome_path("{accession}.fna.gz"),
        processed=outpath("preprocessed", "{accession}.fa")
    output:
        outpath("preprocessing", "{accession}.summary.tsv")
    params:
        accession="{accession}"
    script:
        "../scripts/summarize_preprocessing.py"


rule preprocess_summary:
    input:
        expand(outpath("preprocessing", "{accession}.summary.tsv"), accession=SAMPLE_IDS)
    output:
        outpath("preprocessing", "preprocessing_summary.tsv")
    script:
        "../scripts/merge_tables.py"
