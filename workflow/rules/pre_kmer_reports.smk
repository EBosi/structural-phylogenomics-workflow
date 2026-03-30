rule pre_kmer_summary:
    input:
        assemblies="results/metadata/assemblies.tsv",
        qc="results/qc/qc_summary.tsv",
        preprocessing="results/preprocessing/preprocessing_summary.tsv",
        organelle="results/organelle/organelle_summary.tsv",
        repeats="results/repeats/repeat_annotation_summary.tsv"
    output:
        "results/reports/pre_kmer_summary.tsv"
    script:
        "../scripts/build_pre_kmer_summary.py"


rule pre_kmer_report:
    input:
        assemblies="results/metadata/assemblies.tsv",
        summary="results/reports/pre_kmer_summary.tsv"
    output:
        "results/reports/pre_kmer_report.md"
    params:
        repeat_backend=config["repeat_annotation"]["backend"]
    script:
        "../scripts/render_pre_kmer_report.py"
