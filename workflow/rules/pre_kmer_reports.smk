rule pre_kmer_summary:
    input:
        assemblies=outpath("metadata", "assemblies.tsv"),
        qc=outpath("qc", "qc_summary.tsv"),
        preprocessing=outpath("preprocessing", "preprocessing_summary.tsv"),
        organelle=outpath("organelle", "organelle_summary.tsv"),
        repeats=outpath("repeats", "repeat_annotation_summary.tsv")
    output:
        outpath("reports", "pre_kmer_summary.tsv")
    script:
        "../scripts/build_pre_kmer_summary.py"


rule pre_kmer_report:
    input:
        assemblies=outpath("metadata", "assemblies.tsv"),
        summary=outpath("reports", "pre_kmer_summary.tsv")
    output:
        outpath("reports", "pre_kmer_report.md")
    params:
        repeat_backend=config["repeat_annotation"]["backend"]
    script:
        "../scripts/render_pre_kmer_report.py"
