rule infer_tree:
    input:
        outpath("distances", "{dataset}", "k{k}", "{metric}.tsv")
    output:
        outpath("trees", "{dataset}", "k{k}", "{metric}", "{method}.nwk")
    params:
        method="{method}"
    script:
        "../scripts/infer_tree.py"
