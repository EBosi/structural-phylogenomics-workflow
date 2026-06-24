rule tree_manifest:
    input:
        expand(
            outpath("trees", "{dataset}", "k{k}", "{metric}", "{method}.nwk"),
            dataset=KMER_DATASETS,
            k=K_VALUES,
            metric=DISTANCE_METRICS,
            method=TREE_METHODS,
        )
    output:
        outpath("reports", "tree_manifest.tsv")
    script:
        "../scripts/build_tree_manifest.py"


rule tree_comparisons:
    input:
        manifest=outpath("reports", "tree_manifest.tsv")
    output:
        outpath("reports", "tree_comparisons.tsv")
    script:
        "../scripts/compare_trees.py"
