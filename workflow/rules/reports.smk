rule tree_manifest:
    input:
        expand(
            "results/trees/{dataset}/k{k}/{metric}/{method}.nwk",
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
            method=TREE_METHODS,
        )
    output:
        "results/reports/tree_manifest.tsv"
    script:
        "../scripts/build_tree_manifest.py"


rule tree_comparisons:
    input:
        manifest="results/reports/tree_manifest.tsv"
    output:
        "results/reports/tree_comparisons.tsv"
    script:
        "../scripts/compare_trees.py"
