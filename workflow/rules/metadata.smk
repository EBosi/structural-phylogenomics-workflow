rule fetch_ncbi_assembly_summary_genbank:
    output:
        "resources/ncbi/assembly_summary_genbank.txt"
    params:
        url=config["metadata"]["ncbi_genbank_summary_url"],
        curl=config["downloads"]["curl_path"]
    shell:
        "{params.curl} -fsSL {params.url} -o {output}"


rule fetch_ncbi_assembly_summary_refseq:
    output:
        "resources/ncbi/assembly_summary_refseq.txt"
    params:
        url=config["metadata"]["ncbi_refseq_summary_url"],
        curl=config["downloads"]["curl_path"]
    shell:
        "{params.curl} -fsSL {params.url} -o {output}"


rule resolve_accession_metadata:
    input:
        accession_file=config["metadata"]["accession_file"],
        genbank="resources/ncbi/assembly_summary_genbank.txt",
        refseq="resources/ncbi/assembly_summary_refseq.txt"
    output:
        assemblies="results/metadata/assemblies.tsv",
        organisms="results/metadata/organisms.tsv",
        manifest="results/metadata/download_manifest.tsv"
    params:
        accessions=ACCESSIONS
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
