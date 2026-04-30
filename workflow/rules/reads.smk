rule build_reads_manifest:
    input:
        metadata=config["metadata"]["local_reads_file"]
    output:
        "results/reads/metadata/samples.tsv"
    params:
        required_columns=config["reads"]["metadata_required_columns"],
        fastq_extensions=config["reads"]["fastq_extensions"]
    script:
        "../scripts/build_reads_manifest.py"


rule read_qc:
    input:
        manifest="results/reads/metadata/samples.tsv"
    output:
        "results/reads/qc/{sample}.tsv"
    script:
        "../scripts/compute_read_qc.py"


rule summarize_read_qc:
    input:
        manifest="results/reads/metadata/samples.tsv",
        sample_reports=expand("results/reads/qc/{sample}.tsv", sample=READ_SAMPLE_IDS)
    output:
        "results/reads/qc/summary.tsv"
    script:
        "../scripts/summarize_read_qc.py"


rule model_read_k_response:
    input:
        manifest="results/reads/metadata/samples.tsv",
        qc="results/reads/qc/{sample}.tsv",
        histogram="results/reads/histograms/{sample}.tsv"
    output:
        "results/reads/models/{sample}.tsv"
    params:
        candidate_k_values=config["reads"]["candidate_k_values"],
        min_supported_k=config["reads"]["min_supported_k"],
        low_coverage_threshold=config["reads"]["low_coverage_threshold"],
        medium_coverage_threshold=config["reads"]["medium_coverage_threshold"],
        high_coverage_threshold=config["reads"]["high_coverage_threshold"],
        min_mean_read_length=config["reads"]["min_mean_read_length"],
        include_with_caution_coverage_threshold=config["reads"]["include_with_caution_coverage_threshold"],
        min_retained_fraction=config["reads"]["min_retained_fraction"],
        max_singleton_fraction=config["reads"]["max_singleton_fraction"],
        max_high_copy_fraction=config["reads"]["max_high_copy_fraction"],
        min_mean_quality=config["reads"]["min_mean_quality"],
        max_n_fraction=config["reads"]["max_n_fraction"],
        exclude_pair_mismatch=config["reads"]["exclude_pair_mismatch"]
    script:
        "../scripts/model_read_k_response.py"


rule compute_read_kmer_histogram:
    input:
        manifest="results/reads/metadata/samples.tsv",
        qc="results/reads/qc/{sample}.tsv"
    output:
        "results/reads/histograms/{sample}.tsv"
    params:
        candidate_k_values=config["reads"]["candidate_k_values"],
        max_reads_per_file=config["reads"]["histogram_max_reads_per_file"],
        max_count_bin=config["reads"]["histogram_max_count_bin"]
    script:
        "../scripts/compute_read_kmer_histogram.py"


rule summarize_read_histograms:
    input:
        manifest="results/reads/metadata/samples.tsv",
        histogram_reports=expand("results/reads/histograms/{sample}.tsv", sample=READ_SAMPLE_IDS)
    output:
        "results/reads/histograms/summary.tsv"
    script:
        "../scripts/summarize_read_histograms.py"


rule summarize_read_models:
    input:
        manifest="results/reads/metadata/samples.tsv",
        model_reports=expand("results/reads/models/{sample}.tsv", sample=READ_SAMPLE_IDS)
    output:
        "results/reads/models/summary.tsv"
    script:
        "../scripts/summarize_read_models.py"


rule score_read_samples:
    input:
        manifest="results/reads/metadata/samples.tsv",
        model_summary="results/reads/models/summary.tsv"
    output:
        "results/reads/sample_scores.tsv"
    script:
        "../scripts/score_read_samples.py"


rule select_dataset_k:
    input:
        sample_scores="results/reads/sample_scores.tsv"
    output:
        "results/reads/dataset_k_selection.tsv"
    params:
        preferred_k_values=config["reads"]["preferred_k_values"],
        min_supported_k=config["reads"]["min_supported_k"],
        support_fraction=config["reads"]["dataset_k_support_fraction"]
    script:
        "../scripts/select_dataset_k.py"


rule infer_read_abundance_band:
    input:
        histogram="results/reads/histograms/{sample}.tsv",
        model="results/reads/models/{sample}.tsv",
        dataset_k="results/reads/dataset_k_selection.tsv"
    output:
        "results/reads/bands/{sample}.tsv"
    script:
        "../scripts/infer_read_abundance_band.py"


rule compute_filtered_read_signature:
    input:
        manifest="results/reads/metadata/samples.tsv",
        band="results/reads/bands/{sample}.tsv",
        dataset_k="results/reads/dataset_k_selection.tsv"
    output:
        "results/reads/sketch/signatures/{sample}.tsv"
    params:
        num_hashes=config["reads"]["sketch_num_hashes"],
        max_reads_per_file=config["reads"]["sketch_max_reads_per_file"],
        exclude_low_confidence_bands=config["reads"]["exclude_low_confidence_bands"]
    script:
        "../scripts/compute_filtered_read_signature.py"


rule compute_read_sketch_distance_matrix:
    input:
        expand("results/reads/sketch/signatures/{sample}.tsv", sample=READ_SAMPLE_IDS)
    output:
        "results/reads/sketch/distances/minhash_jaccard.tsv"
    params:
        samples=READ_SAMPLE_IDS,
        exclude_signature_statuses=[
            "excluded_by_qc",
            "unsupported_dataset_k",
            "low_confidence_band",
            "no_retained_kmers",
        ]
    script:
        "../scripts/compute_sketch_distance_matrix.py"


rule infer_read_sketch_tree:
    input:
        "results/reads/sketch/distances/minhash_jaccard.tsv"
    output:
        "results/reads/sketch/tree.nwk"
    params:
        method="nj"
    script:
        "../scripts/infer_tree_safe.py"


rule render_reads_report:
    input:
        manifest="results/reads/metadata/samples.tsv",
        qc_summary="results/reads/qc/summary.tsv",
        histogram_summary="results/reads/histograms/summary.tsv",
        model_summary="results/reads/models/summary.tsv",
        band_reports=expand("results/reads/bands/{sample}.tsv", sample=READ_SAMPLE_IDS),
        signature_reports=expand("results/reads/sketch/signatures/{sample}.tsv", sample=READ_SAMPLE_IDS),
        sample_scores="results/reads/sample_scores.tsv",
        dataset_k="results/reads/dataset_k_selection.tsv",
        sketch_distance="results/reads/sketch/distances/minhash_jaccard.tsv",
        sketch_tree="results/reads/sketch/tree.nwk"
    output:
        "results/reads/report.md"
    script:
        "../scripts/render_reads_report.py"
