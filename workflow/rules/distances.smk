rule distance_matrix:
    input:
        "results/kmers/matrices/{dataset}/k{k}.tsv"
    output:
        "results/distances/{dataset}/k{k}/{metric}.tsv"
    params:
        metric="{metric}"
    script:
        "../scripts/compute_distance_matrix.py"
