import csv
import json
import urllib.parse
import urllib.request
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


def localize_ftp_path(ftp_path):
    if ftp_path.startswith("ftp://"):
        return "https://" + ftp_path[len("ftp://") :]
    return ftp_path


def record_from_summary_row(row, source_db):
    ftp_path = row.get("ftp_path", "")
    if ftp_path == "na":
        ftp_path = ""
    ftp_path = localize_ftp_path(ftp_path)
    return {
        "datatype": "genome",
        "accession": row["# assembly_accession"],
        "resolved_accession": row["# assembly_accession"],
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
        "local_path": f"data/genomes/{row['# assembly_accession']}.fna.gz",
    }


def esearch_uid(accession, eutils_base, timeout):
    term = urllib.parse.quote(f"{accession}[Assembly Accession]")
    url = f"{eutils_base}/esearch.fcgi?db=assembly&term={term}&retmode=json"
    payload = json.loads(urllib.request.urlopen(url, timeout=timeout).read().decode("utf-8"))
    idlist = payload["esearchresult"]["idlist"]
    if not idlist:
        raise ValueError(f"NCBI ESearch returned no assembly UID for accession {accession}")
    return idlist[0]


def fetch_esummary_records(accessions, eutils_base, timeout):
    accession_to_uid = {accession: esearch_uid(accession, eutils_base, timeout) for accession in accessions}
    uids = [accession_to_uid[accession] for accession in accessions]
    summary_url = f"{eutils_base}/esummary.fcgi?db=assembly&id={','.join(uids)}&retmode=json"
    payload = json.loads(urllib.request.urlopen(summary_url, timeout=timeout).read().decode("utf-8"))

    records = {}
    for accession in accessions:
        uid = accession_to_uid[accession]
        raw = payload["result"][uid]
        ftp_refseq = localize_ftp_path(raw.get("ftppath_refseq", "") or "")
        ftp_genbank = localize_ftp_path(raw.get("ftppath_genbank", "") or "")
        ftp_path = ftp_refseq or ftp_genbank
        source_db = "refseq" if ftp_refseq else "genbank"
        records[accession] = {
            "datatype": "genome",
            "accession": accession,
            "resolved_accession": raw.get("assemblyaccession", accession),
            "bioproject": raw.get("bioprojectaccn", "") or raw.get("bioprojectid", ""),
            "biosample": raw.get("biosampleaccn", ""),
            "organism_name": raw.get("organism", ""),
            "species_taxid": str(raw.get("speciestaxid", "") or ""),
            "taxid": str(raw.get("taxid", "") or ""),
            "infraspecific_name": raw.get("infraspecieslist", "") or "",
            "isolate": raw.get("isolate", ""),
            "assembly_name": raw.get("assemblyname", ""),
            "assembly_level": raw.get("assemblystatus", ""),
            "version_status": raw.get("versionstatus", ""),
            "release_type": raw.get("releasetype", ""),
            "genome_rep": raw.get("genomerep", ""),
            "seq_rel_date": raw.get("seqreldate", ""),
            "submitter": raw.get("submitterorganization", ""),
            "refseq_category": raw.get("refseq_category", "") or raw.get("refseqcategory", ""),
            "source_db": source_db,
            "ftp_path": ftp_path,
            "ftp_url": build_ftp_url(ftp_path),
            "local_path": f"data/genomes/{accession}.fna.gz",
        }
    return records


def load_local_genome_records(path):
    records = []
    local_path = Path(path)
    if not local_path.exists():
        return records

    with local_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return records
        required = {"accession", "local_path"}
        missing = required - set(reader.fieldnames)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise ValueError(
                f"Local genome metadata file {local_path} is missing required columns: {missing_str}"
            )
        for row in reader:
            accession = row.get("accession", "").strip()
            fasta_path = row.get("local_path", "").strip()
            if not accession:
                continue
            if not fasta_path:
                raise ValueError(f"Local genome entry {accession} is missing a local_path value")
            records.append(
                {
                    "datatype": row.get("datatype", "genome").strip() or "genome",
                    "accession": accession,
                    "resolved_accession": row.get("resolved_accession", accession).strip() or accession,
                    "bioproject": row.get("bioproject", "").strip(),
                    "biosample": row.get("biosample", "").strip(),
                    "organism_name": row.get("organism_name", accession).strip() or accession,
                    "species_taxid": row.get("species_taxid", "").strip(),
                    "taxid": row.get("taxid", "").strip(),
                    "infraspecific_name": row.get("infraspecific_name", "").strip(),
                    "isolate": row.get("isolate", "").strip(),
                    "assembly_name": row.get("assembly_name", "").strip(),
                    "assembly_level": row.get("assembly_level", "").strip(),
                    "version_status": row.get("version_status", "local").strip() or "local",
                    "release_type": row.get("release_type", "local").strip() or "local",
                    "genome_rep": row.get("genome_rep", "full").strip() or "full",
                    "seq_rel_date": row.get("seq_rel_date", "").strip(),
                    "submitter": row.get("submitter", "").strip(),
                    "refseq_category": row.get("refseq_category", "").strip(),
                    "source_db": row.get("source_db", "local").strip() or "local",
                    "ftp_path": "",
                    "ftp_url": "",
                    "local_path": fasta_path,
                }
            )
    return records


requested = list(snakemake.params.accessions)
timeout = int(getattr(snakemake.params, "request_timeout", 60))
eutils_base = getattr(snakemake.params, "eutils_base", "")
local_records = load_local_genome_records(snakemake.input.local_genomes)

if requested and hasattr(snakemake.input, "genbank") and hasattr(snakemake.input, "refseq"):
    records = {}
    records.update(load_summary_table(snakemake.input.genbank, "genbank"))
    records.update(load_summary_table(snakemake.input.refseq, "refseq"))
elif requested:
    if not eutils_base:
        raise ValueError("eutils_base must be provided when resolving accessions from NCBI E-utilities")
    records = fetch_esummary_records(requested, eutils_base, timeout)
else:
    records = {}

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

assemblies.extend(local_records)

assembly_output = Path(snakemake.output.assemblies)
assembly_output.parent.mkdir(parents=True, exist_ok=True)

assembly_fields = [
        "datatype",
        "accession",
        "resolved_accession",
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
