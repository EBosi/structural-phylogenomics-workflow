rule distance_matrix:
    input:
        outpath("kmers", "matrices", "{dataset}", "k{k}.tsv")
    output:
        outpath("distances", "{dataset}", "k{k}", "{metric}.tsv")
    params:
        metric="{metric}"
    script:
        "../scripts/compute_distance_matrix.py"
