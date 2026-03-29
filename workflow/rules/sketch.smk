def sketch_input_for_dataset(wildcards):
    if wildcards.dataset == "unmasked":
        return f"results/organelle/filtered/{wildcards.accession}.fa"
    if wildcards.dataset == "masked":
        return f"results/repeats/masked/{wildcards.accession}.fa"
    raise ValueError(f"Unsupported dataset: {wildcards.dataset}")


rule sketch_signature_per_sample:
    input:
        sketch_input_for_dataset
    output:
        "results/sketch/signatures/{dataset}/k{k}/{accession}.tsv"
    params:
        accession="{accession}",
        dataset="{dataset}",
        k="{k}",
        num_hashes=config["sketch"]["num_hashes"],
        canonical=config["sketch"]["canonical"]
    threads: 1
    script:
        "../scripts/compute_minhash_signature.py"


rule sketch_distance_matrix:
    input:
        lambda wc: expand(
            "results/sketch/signatures/{dataset}/k{k}/{accession}.tsv",
            dataset=wc.dataset,
            k=wc.k,
            accession=SAMPLE_IDS,
        )
    output:
        "results/sketch/distances/{dataset}/k{k}/minhash_jaccard.tsv"
    params:
        accessions=SAMPLE_IDS
    threads: 1
    script:
        "../scripts/compute_sketch_distance_matrix.py"


rule infer_sketch_tree:
    input:
        "results/sketch/distances/{dataset}/k{k}/minhash_jaccard.tsv"
    output:
        "results/sketch/trees/{dataset}/k{k}/{method}.nwk"
    params:
        method="{method}"
    threads: 1
    script:
        "../scripts/infer_tree.py"
