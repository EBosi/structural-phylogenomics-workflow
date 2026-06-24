def kmer_input_for_dataset(wildcards):
    if wildcards.dataset == "unmasked":
        return outpath("organelle", "filtered", f"{wildcards.accession}.fa")
    if wildcards.dataset == "masked":
        return outpath("repeats", "masked", f"{wildcards.accession}.fa")
    raise ValueError(f"Unsupported dataset: {wildcards.dataset}")


rule kmer_spectrum_per_sample:
    input:
        kmer_input_for_dataset
    output:
        outpath("kmers", "spectra", "{dataset}", "k{k}", "{accession}.tsv")
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
            outpath("kmers", "spectra", "{dataset}", "k{k}", "{accession}.tsv"),
            dataset=wc.dataset,
            k=wc.k,
            accession=SAMPLE_IDS,
        )
    output:
        outpath("kmers", "matrices", "{dataset}", "k{k}.tsv")
    params:
        accessions=SAMPLE_IDS,
        value_column=config["kmers"]["normalization"]
    script:
        "../scripts/build_kmer_matrix.py"
