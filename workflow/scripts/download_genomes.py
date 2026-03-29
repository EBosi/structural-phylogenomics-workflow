import csv
import subprocess
from pathlib import Path


manifest_path = Path(snakemake.input.manifest)
curl_path = snakemake.params.curl

rows = []
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows.extend(reader)

for row in rows:
    ftp_url = row["ftp_url"].strip()
    local_path = Path(row["local_path"])
    accession = row["accession"]
    local_path.parent.mkdir(parents=True, exist_ok=True)

    if not ftp_url:
        if local_path.exists():
            continue
        raise ValueError(
            f"No ftp_url available for accession {accession} and local file {local_path} does not exist"
        )

    subprocess.run(
        [curl_path, "-fsSL", ftp_url, "-o", str(local_path)],
        check=True,
    )
