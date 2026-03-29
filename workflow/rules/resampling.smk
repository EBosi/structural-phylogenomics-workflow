def resampling_input_for_dataset(wildcards):
    if wildcards.dataset == "unmasked":
        return f"results/organelle/filtered/{wildcards.accession}.fa"
    if wildcards.dataset == "masked":
        return f"results/repeats/masked/{wildcards.accession}.fa"
    raise ValueError(f"Unsupported dataset: {wildcards.dataset}")


rule unit_kmer_spectra_windows:
    input:
        resampling_input_for_dataset
    output:
        matrix="results/resampling/units/windows/{dataset}/k{k}/{accession}.npy",
        meta="results/resampling/units/windows/{dataset}/k{k}/{accession}.meta.tsv"
    params:
        accession="{accession}",
        dataset="{dataset}",
        k="{k}",
        unit_type="window",
        window_size=config["resampling"]["window_size"],
        step_size=config["resampling"]["step_size"],
        min_window_length=config["resampling"]["min_window_length"],
        canonical=config["kmers"]["canonical"]
    threads: 1
    resources:
        mem_mb=4096
    script:
        "../scripts/compute_partitioned_kmer_spectra.py"


rule unit_kmer_spectra_contigs:
    input:
        resampling_input_for_dataset
    output:
        matrix="results/resampling/units/contigs/{dataset}/k{k}/{accession}.npy",
        meta="results/resampling/units/contigs/{dataset}/k{k}/{accession}.meta.tsv"
    params:
        accession="{accession}",
        dataset="{dataset}",
        k="{k}",
        unit_type="contig",
        window_size=0,
        step_size=0,
        min_window_length=1,
        canonical=config["kmers"]["canonical"]
    threads: 1
    resources:
        mem_mb=8192
    script:
        "../scripts/compute_partitioned_kmer_spectra.py"


rule bootstrap_support:
    input:
        tree="results/trees/{dataset}/k{k}/{metric}/{method}.nwk",
        matrices=lambda wc: expand(
            "results/resampling/units/windows/{dataset}/k{k}/{accession}.npy",
            dataset=wc.dataset,
            k=wc.k,
            accession=ACCESSIONS,
        ),
        metadata=lambda wc: expand(
            "results/resampling/units/windows/{dataset}/k{k}/{accession}.meta.tsv",
            dataset=wc.dataset,
            k=wc.k,
            accession=ACCESSIONS,
        )
    output:
        support="results/resampling/bootstrap/{dataset}/k{k}/{metric}/{method}.support.tsv",
        summary="results/resampling/bootstrap/{dataset}/k{k}/{metric}/{method}.summary.tsv"
    params:
        mode="bootstrap",
        metric="{metric}",
        method="{method}",
        replicates=config["resampling"]["bootstrap_replicates"],
        fraction=1.0,
        seed=config["resampling"]["random_seed"]
    threads: 1
    resources:
        mem_mb=8192
    script:
        "../scripts/resample_tree_support.py"


rule jackknife_support:
    input:
        tree="results/trees/{dataset}/k{k}/{metric}/{method}.nwk",
        matrices=lambda wc: expand(
            "results/resampling/units/contigs/{dataset}/k{k}/{accession}.npy",
            dataset=wc.dataset,
            k=wc.k,
            accession=ACCESSIONS,
        ),
        metadata=lambda wc: expand(
            "results/resampling/units/contigs/{dataset}/k{k}/{accession}.meta.tsv",
            dataset=wc.dataset,
            k=wc.k,
            accession=ACCESSIONS,
        )
    output:
        support="results/resampling/jackknife/{dataset}/k{k}/{metric}/{method}.support.tsv",
        summary="results/resampling/jackknife/{dataset}/k{k}/{metric}/{method}.summary.tsv"
    params:
        mode="jackknife",
        metric="{metric}",
        method="{method}",
        replicates=config["resampling"]["jackknife_replicates"],
        fraction=config["resampling"]["jackknife_fraction"],
        seed=config["resampling"]["random_seed"]
    threads: 1
    resources:
        mem_mb=8192
    script:
        "../scripts/resample_tree_support.py"
