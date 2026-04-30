import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from repeat_utils import (
    combine_record_intervals,
    interval_map_to_records,
    load_interval_text,
    parse_repeatmasker_out_records,
    write_interval_details,
    write_interval_text,
)
from io_utils import read_fasta


def backend_components(backend):
    components = {component.strip().lower() for component in backend.split("+") if component.strip()}
    supported = {"dustmasker", "repeatmasker", "repeatmodeler"}
    unsupported = components - supported
    if unsupported:
        supported_str = ", ".join(sorted(supported))
        unsupported_str = ", ".join(sorted(unsupported))
        raise ValueError(f"Unsupported repeat backend component(s): {unsupported_str}. Supported: {supported_str}")
    if not components:
        raise ValueError("repeat_annotation.backend must not be empty")
    return components


def resolve_executable(path_value, executable_name):
    if path_value:
        configured = Path(path_value)
        if configured.is_absolute():
            if not configured.exists():
                raise FileNotFoundError(
                    f"Configured executable for {executable_name} does not exist: {configured}"
                )
            return str(configured)
        found = shutil.which(path_value)
        if found:
            return found
    found = shutil.which(executable_name)
    if found:
        return found
    raise FileNotFoundError(
        f"Could not resolve executable for {executable_name}. Set an absolute path in config/config.yaml."
    )


def repeatmasker_env(repeatmasker_path):
    env = os.environ.copy()
    repeatmasker_bin = str(Path(repeatmasker_path).resolve().parent)
    env["PATH"] = f"{repeatmasker_bin}:{env.get('PATH', '')}"
    return env


def run_dustmasker(input_fasta, output_path, params):
    dustmasker = resolve_executable(params.dustmasker, "dustmasker")
    subprocess.run(
        [
            dustmasker,
            "-in",
            input_fasta,
            "-out",
            output_path,
            "-outfmt",
            "interval",
            "-window",
            str(params.window),
            "-level",
            str(params.level),
            "-linker",
            str(params.linker),
        ],
        check=True,
    )
    intervals = load_interval_text(output_path)
    return interval_map_to_records(intervals, "dustmasker", "dustmasker", "Low_complexity")


def copy_repeatmasker_reports(workdir, raw_dir, source):
    for suffix in (".out", ".tbl"):
        for path in sorted(Path(workdir).glob(f"*{suffix}")):
            shutil.copy2(path, raw_dir / f"{source}{suffix}")
            break


def run_repeatmasker(input_fasta, workdir, params, source, species="", library=""):
    repeatmasker = resolve_executable(params.repeatmasker, "RepeatMasker")
    repeatmasker_fasta, id_map = write_repeatmasker_safe_fasta(input_fasta, workdir / "repeatmasker_input.fa")
    command = [repeatmasker, "-dir", str(workdir), "-pa", str(params.repeatmasker_threads)]
    if params.repeatmasker_engine:
        command.extend(["-engine", params.repeatmasker_engine])
    if species:
        command.extend(["-species", species])
    if library:
        command.extend(["-lib", library])
    if params.repeatmasker_extra_args:
        command.extend(shlex.split(params.repeatmasker_extra_args))
    command.append(str(repeatmasker_fasta))

    subprocess.run(command, check=True, env=repeatmasker_env(repeatmasker))

    out_files = sorted(Path(workdir).glob("*.out"))
    if not out_files:
        raise FileNotFoundError("RepeatMasker completed but no .out file was produced")
    records = parse_repeatmasker_out_records(out_files[0], source=source)
    for record in records:
        record["seq_id"] = id_map.get(record["seq_id"], record["seq_id"])
    return records


def write_repeatmasker_safe_fasta(input_fasta, output_fasta):
    id_map = {}
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w") as handle:
        for index, (header, sequence) in enumerate(read_fasta(input_fasta), start=1):
            original_id = header.split()[0]
            safe_id = f"seq{index:08d}"
            id_map[safe_id] = original_id
            handle.write(f">{safe_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    return output_fasta, id_map


def write_fasta_prefix(input_fasta, output_fasta, max_bases):
    written = 0
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w") as handle:
        for header, sequence in read_fasta(input_fasta):
            if written >= max_bases:
                break
            remaining = max_bases - written
            trimmed = sequence[:remaining]
            if not trimmed:
                break
            handle.write(f">{header}\n")
            for start in range(0, len(trimmed), 80):
                handle.write(trimmed[start : start + 80] + "\n")
            written += len(trimmed)
    if written == 0:
        raise ValueError(f"No sequence bases available for bounded repeat annotation: {input_fasta}")
    return output_fasta


def denovo_library_input(snakemake_input):
    value = getattr(snakemake_input, "denovo_library", [])
    if isinstance(value, str):
        return value
    if value:
        return value[0]
    return ""


backend = snakemake.params.backend
components = backend_components(backend)

input_fasta = str(snakemake.input.fasta)
output_intervals = Path(snakemake.output.intervals)
output_details = Path(snakemake.output.details)
raw_dir = Path(snakemake.output.raw_dir)
output_intervals.parent.mkdir(parents=True, exist_ok=True)
output_details.parent.mkdir(parents=True, exist_ok=True)
if raw_dir.exists():
    shutil.rmtree(raw_dir)
raw_dir.mkdir(parents=True, exist_ok=True)

denovo_library = denovo_library_input(snakemake.input)
needs_known_repeatmasker = "repeatmasker" in components
needs_denovo_repeatmasker = "repeatmodeler" in components
has_known_selector = bool(snakemake.params.repeatmasker_species or snakemake.params.repeatmasker_library)

if needs_known_repeatmasker and not has_known_selector and bool(snakemake.params.repeatmasker_require_species_or_library):
    raise ValueError(
        "RepeatMasker backend requested without repeatmasker_species or repeatmasker_library. "
        "Do not rely on RepeatMasker's default database selection for non-model eukaryotes."
    )
if needs_denovo_repeatmasker and not denovo_library:
    raise ValueError("RepeatModeler backend requested but no de novo library input was provided")

records = []

with tempfile.TemporaryDirectory(prefix="repeat_annotation_") as temp_dir:
    temp_path = Path(temp_dir)
    max_input_bases = int(getattr(snakemake.params, "max_input_bases", 0) or 0)
    annotation_fasta = input_fasta
    if max_input_bases > 0:
        annotation_fasta = str(write_fasta_prefix(input_fasta, temp_path / "bounded_input.fa", max_input_bases))

    if "dustmasker" in components:
        dust_intervals = temp_path / "dustmasker.intervals.txt"
        dust_records = run_dustmasker(annotation_fasta, str(dust_intervals), snakemake.params)
        shutil.copy2(dust_intervals, raw_dir / "dustmasker.intervals.txt")
        records.extend(dust_records)

    if needs_known_repeatmasker:
        known_dir = temp_path / "repeatmasker_known"
        known_dir.mkdir()
        source = "repeatmasker_custom" if snakemake.params.repeatmasker_library else "repeatmasker_known"
        known_records = run_repeatmasker(
            annotation_fasta,
            known_dir,
            snakemake.params,
            source=source,
            species=snakemake.params.repeatmasker_species,
            library=snakemake.params.repeatmasker_library,
        )
        copy_repeatmasker_reports(known_dir, raw_dir, source)
        records.extend(known_records)

    if needs_denovo_repeatmasker:
        denovo_dir = temp_path / "repeatmasker_denovo"
        denovo_dir.mkdir()
        denovo_records = run_repeatmasker(
            annotation_fasta,
            denovo_dir,
            snakemake.params,
            source="repeatmasker_denovo",
            library=denovo_library,
        )
        copy_repeatmasker_reports(denovo_dir, raw_dir, "repeatmasker_denovo")
        records.extend(denovo_records)

combined = combine_record_intervals(records)
write_interval_text(combined, output_intervals)
write_interval_details(records, output_details)
