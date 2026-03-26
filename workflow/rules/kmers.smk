def kmer_input_for_dataset(wildcards):
    if wildcards.dataset == "unmasked":
        return f"results/preprocessed/{wildcards.accession}.fa"
    if wildcards.dataset == "masked":
        return f"results/repeats/masked/{wildcards.accession}.fa"
    raise ValueError(f"Unsupported dataset: {wildcards.dataset}")


rule kmer_spectrum_per_sample:
    input:
        kmer_input_for_dataset
    output:
        "results/kmers/spectra/{dataset}/k{k}/{accession}.tsv"
    params:
        accession="{accession}",
        dataset="{dataset}",
        k="{k}",
        canonical=config["kmers"]["canonical"],
        normalization=config["kmers"]["normalization"]
    script:
        "../scripts/compute_kmer_spectrum.py"


rule kmer_feature_matrix:
    input:
        lambda wc: expand(
            "results/kmers/spectra/{dataset}/k{k}/{accession}.tsv",
            dataset=wc.dataset,
            k=wc.k,
            accession=ACCESSIONS,
        )
    output:
        "results/kmers/matrices/{dataset}/k{k}.tsv"
    params:
        accessions=ACCESSIONS,
        value_column=config["kmers"]["normalization"]
    script:
        "../scripts/build_kmer_matrix.py"
