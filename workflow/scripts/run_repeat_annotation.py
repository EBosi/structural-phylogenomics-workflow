import shlex
import subprocess
import tempfile
import os
from pathlib import Path

from repeat_utils import combine_interval_sets, load_interval_text, parse_repeatmasker_out, write_interval_text


def run_dustmasker(input_fasta, output_path, params):
    subprocess.run(
        [
            params.dustmasker,
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
    return load_interval_text(output_path)


def run_repeatmasker(input_fasta, workdir, params):
    command = [params.repeatmasker, "-dir", str(workdir), "-pa", str(params.repeatmasker_threads)]
    if params.repeatmasker_engine:
        command.extend(["-engine", params.repeatmasker_engine])
    if params.repeatmasker_species:
        command.extend(["-species", params.repeatmasker_species])
    if params.repeatmasker_library:
        command.extend(["-lib", params.repeatmasker_library])
    if params.repeatmasker_extra_args:
        command.extend(shlex.split(params.repeatmasker_extra_args))
    command.append(input_fasta)

    env = os.environ.copy()
    repeatmasker_bin = str(Path(params.repeatmasker).resolve().parent)
    env["PATH"] = f"{repeatmasker_bin}:{env.get('PATH', '')}"

    subprocess.run(command, check=True, env=env)

    out_files = sorted(Path(workdir).glob("*.out"))
    if not out_files:
        raise FileNotFoundError("RepeatMasker completed but no .out file was produced")
    return parse_repeatmasker_out(out_files[0])


backend = snakemake.params.backend
supported_backends = {"dustmasker", "repeatmasker", "dustmasker+repeatmasker"}
if backend not in supported_backends:
    supported_str = ", ".join(sorted(supported_backends))
    raise ValueError(f"Unsupported repeat backend '{backend}'. Supported values: {supported_str}")

input_fasta = str(snakemake.input[0])
output_intervals = Path(snakemake.output[0])
output_intervals.parent.mkdir(parents=True, exist_ok=True)

interval_sets = []

with tempfile.TemporaryDirectory(prefix="repeat_annotation_") as temp_dir:
    temp_path = Path(temp_dir)

    if backend in {"dustmasker", "dustmasker+repeatmasker"}:
        dust_intervals = temp_path / "dustmasker.intervals.txt"
        interval_sets.append(run_dustmasker(input_fasta, str(dust_intervals), snakemake.params))

    if backend in {"repeatmasker", "dustmasker+repeatmasker"}:
        interval_sets.append(run_repeatmasker(input_fasta, temp_path, snakemake.params))

combined = combine_interval_sets(interval_sets)
write_interval_text(combined, output_intervals)
