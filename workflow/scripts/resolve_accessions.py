import csv
from collections import defaultdict
from pathlib import Path


def load_summary_table(path, source_db):
    records = {}
    header = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("# assembly_accession"):
                header = line.rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            if header is None:
                raise ValueError(f"Header not found in {path}")
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            accession = row["# assembly_accession"]
            ftp_path = row.get("ftp_path", "")
            if ftp_path == "na":
                ftp_path = ""
            records[accession] = {
                "datatype": "genome",
                "accession": accession,
                "bioproject": row.get("bioproject", ""),
                "biosample": row.get("biosample", ""),
                "organism_name": row.get("organism_name", ""),
                "species_taxid": row.get("species_taxid", ""),
                "taxid": row.get("taxid", ""),
                "infraspecific_name": row.get("infraspecific_name", ""),
                "isolate": row.get("isolate", ""),
                "assembly_name": row.get("asm_name", ""),
                "assembly_level": row.get("assembly_level", ""),
                "version_status": row.get("version_status", ""),
                "release_type": row.get("release_type", ""),
                "genome_rep": row.get("genome_rep", ""),
                "seq_rel_date": row.get("seq_rel_date", ""),
                "submitter": row.get("submitter", ""),
                "refseq_category": row.get("refseq_category", ""),
                "source_db": source_db,
                "ftp_path": ftp_path,
                "ftp_url": build_ftp_url(ftp_path),
                "local_path": f"data/genomes/{accession}.fna.gz",
            }
    return records


def build_ftp_url(ftp_path):
    if not ftp_path:
        return ""
    basename = ftp_path.rstrip("/").split("/")[-1]
    return f"{ftp_path}/{basename}_genomic.fna.gz"


requested = list(snakemake.params.accessions)
records = {}
records.update(load_summary_table(snakemake.input.genbank, "genbank"))
records.update(load_summary_table(snakemake.input.refseq, "refseq"))

assemblies = []
missing = []
for accession in requested:
    record = records.get(accession)
    if record is None:
        missing.append(accession)
        continue
    assemblies.append(record)

if missing:
    missing_str = ", ".join(missing)
    raise ValueError(f"Could not resolve the following accessions from NCBI summary tables: {missing_str}")

assembly_output = Path(snakemake.output.assemblies)
assembly_output.parent.mkdir(parents=True, exist_ok=True)

assembly_fields = [
    "datatype",
    "accession",
    "organism_name",
    "taxid",
    "species_taxid",
    "bioproject",
    "biosample",
    "infraspecific_name",
    "isolate",
    "assembly_name",
    "assembly_level",
    "version_status",
    "release_type",
    "genome_rep",
    "seq_rel_date",
    "submitter",
    "refseq_category",
    "source_db",
    "ftp_path",
    "ftp_url",
    "local_path",
]

with assembly_output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=assembly_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(assemblies)

organism_groups = defaultdict(list)
for row in assemblies:
    organism_groups[(row["organism_name"], row["taxid"], row["species_taxid"])].append(row)

organism_rows = []
for (organism_name, taxid, species_taxid), group in sorted(organism_groups.items()):
    organism_rows.append(
        {
            "organism_name": organism_name,
            "taxid": taxid,
            "species_taxid": species_taxid,
            "n_assemblies": len(group),
            "accessions": ",".join(sorted(item["accession"] for item in group)),
            "assembly_levels": ",".join(sorted({item["assembly_level"] for item in group if item["assembly_level"]})),
            "source_dbs": ",".join(sorted({item["source_db"] for item in group if item["source_db"]})),
        }
    )

with Path(snakemake.output.organisms).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "organism_name",
            "taxid",
            "species_taxid",
            "n_assemblies",
            "accessions",
            "assembly_levels",
            "source_dbs",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(organism_rows)

manifest_rows = [
    {
        "accession": row["accession"],
        "ftp_url": row["ftp_url"],
        "local_path": row["local_path"],
        "source_db": row["source_db"],
    }
    for row in assemblies
]

with Path(snakemake.output.manifest).open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["accession", "ftp_url", "local_path", "source_db"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)
