rule resolve_accession_metadata:
    input:
        accession_file=config["metadata"]["accession_file"]
    output:
        assemblies="results/metadata/assemblies.tsv",
        organisms="results/metadata/organisms.tsv",
        manifest="results/metadata/download_manifest.tsv"
    params:
        accessions=ACCESSIONS,
        eutils_base=config["metadata"]["ncbi_eutils_base"],
        request_timeout=config["metadata"]["request_timeout_seconds"]
    script:
        "../scripts/resolve_accessions.py"


rule download_genomes:
    input:
        manifest="results/metadata/download_manifest.tsv"
    output:
        expand("data/genomes/{accession}.fna.gz", accession=ACCESSIONS)
    params:
        curl=config["downloads"]["curl_path"]
    script:
        "../scripts/download_genomes.py"
