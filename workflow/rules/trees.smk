rule infer_tree:
    input:
        "results/distances/{dataset}/k{k}/{metric}.tsv"
    output:
        "results/trees/{dataset}/k{k}/{metric}/{method}.nwk"
    params:
        method="{method}"
    script:
        "../scripts/infer_tree.py"
