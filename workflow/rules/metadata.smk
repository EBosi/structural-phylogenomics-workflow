rule resolve_accession_metadata:
    input:
        accession_file=config["metadata"]["accession_file"],
        local_genomes=config["metadata"]["local_genomes_file"]
    output:
        assemblies=outpath("metadata", "assemblies.tsv"),
        organisms=outpath("metadata", "organisms.tsv"),
        manifest=outpath("metadata", "download_manifest.tsv")
    params:
        accessions=REMOTE_ACCESSIONS,
        eutils_base=config["metadata"]["ncbi_eutils_base"],
        request_timeout=config["metadata"]["request_timeout_seconds"],
        local_genomes_from_dir=LOCAL_DIR_RECORDS,
        genome_root=GENOME_ROOT
    script:
        "../scripts/resolve_accessions.py"


rule stage_local_genome_from_dir:
    wildcard_constraints:
        accession=LOCAL_DIR_SAMPLE_PATTERN
    input:
        lambda wildcards: LOCAL_DIR_FASTAS[wildcards.accession]
    output:
        genome_path("{accession}.fna.gz")
    script:
        "../scripts/stage_local_genome.py"


rule download_genomes:
    input:
        manifest=outpath("metadata", "download_manifest.tsv")
    output:
        expand(genome_path("{accession}.fna.gz"), accession=REMOTE_ACCESSIONS)
    params:
        curl=config["downloads"]["curl_path"]
    script:
        "../scripts/download_genomes.py"
