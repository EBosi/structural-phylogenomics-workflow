import csv
import re
from pathlib import Path


SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def _normalise_path(path_value: str, metadata_path: Path) -> str:
    if not path_value:
        return ""
    path = Path(path_value)
    if path.is_absolute():
        return str(path)

    cwd_relative = path.resolve()
    metadata_relative = (metadata_path.parent / path).resolve()
    path = cwd_relative if cwd_relative.exists() or not metadata_relative.exists() else metadata_relative
    return str(path)


metadata_path = Path(snakemake.input.metadata)
required_columns = list(snakemake.params.required_columns)
allowed_extensions = tuple(extension.lower() for extension in snakemake.params.fastq_extensions)
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

if not metadata_path.exists():
    raise FileNotFoundError(f"Reads metadata file not found: {metadata_path}")

rows: list[dict[str, str]] = []
with metadata_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"Reads metadata file is empty: {metadata_path}")
    missing = [column for column in required_columns if column not in reader.fieldnames]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Missing required columns in {metadata_path}: {missing_str}")

    seen_sample_ids: set[str] = set()
    for row in reader:
        sample_id = row["sample_id"].strip()
        if not sample_id:
            continue
        if not SAFE_SAMPLE_ID.match(sample_id):
            raise ValueError(
                f"Read sample_id is not safe for filesystem wildcards: {sample_id!r}. "
                "Use only letters, numbers, '.', '_' and '-'."
            )
        if sample_id in seen_sample_ids:
            raise ValueError(f"Duplicated read sample_id in {metadata_path}: {sample_id}")
        seen_sample_ids.add(sample_id)

        read1 = _normalise_path(row["read1"].strip(), metadata_path)
        read2 = _normalise_path(row["read2"].strip(), metadata_path)
        if not read1:
            raise ValueError(f"Sample {sample_id} is missing read1")
        if not read1.lower().endswith(allowed_extensions):
            raise ValueError(f"Sample {sample_id} read1 does not look like FASTQ: {read1}")
        if read2 and not read2.lower().endswith(allowed_extensions):
            raise ValueError(f"Sample {sample_id} read2 does not look like FASTQ: {read2}")
        if not Path(read1).exists():
            raise FileNotFoundError(f"Sample {sample_id} read1 not found: {read1}")
        if read2 and not Path(read2).exists():
            raise FileNotFoundError(f"Sample {sample_id} read2 not found: {read2}")
        genome_size = row["estimated_genome_size_bp"].strip()
        if genome_size:
            try:
                if float(genome_size) <= 0:
                    raise ValueError
            except ValueError as exc:
                raise ValueError(
                    f"Sample {sample_id} estimated_genome_size_bp must be a positive number"
                ) from exc

        rows.append(
            {
                "sample_id": sample_id,
                "species": row["species"].strip(),
                "read1": read1,
                "read2": read2,
                "estimated_genome_size_bp": row["estimated_genome_size_bp"].strip(),
                "platform": row["platform"].strip(),
                "notes": row["notes"].strip(),
            }
        )

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "species",
            "read1",
            "read2",
            "estimated_genome_size_bp",
            "platform",
            "notes",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
