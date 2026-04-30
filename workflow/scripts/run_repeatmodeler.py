import os
import shlex
import shutil
import subprocess
from pathlib import Path


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


input_fasta = str(snakemake.input.fasta)
output_library = Path(snakemake.output.library)
output_summary = Path(snakemake.output.summary)
workdir = Path(snakemake.output.workdir)
if workdir.exists():
    shutil.rmtree(workdir)
workdir.mkdir(parents=True, exist_ok=True)
output_library.parent.mkdir(parents=True, exist_ok=True)

build_database = resolve_executable(snakemake.params.build_database, "BuildDatabase")
repeatmodeler = resolve_executable(snakemake.params.repeatmodeler, "RepeatModeler")
database_name = workdir / snakemake.params.database_name

env = os.environ.copy()
tool_bins = {
    str(Path(build_database).resolve().parent),
    str(Path(repeatmodeler).resolve().parent),
}
env["PATH"] = ":".join(sorted(tool_bins)) + f":{env.get('PATH', '')}"

subprocess.run(
    [build_database, "-name", str(database_name), input_fasta],
    check=True,
    cwd=workdir,
    env=env,
)

command = [repeatmodeler, "-database", str(database_name), "-threads", str(snakemake.threads)]
if bool(snakemake.params.ltr_struct):
    command.append("-LTRStruct")
if snakemake.params.extra_args:
    command.extend(shlex.split(str(snakemake.params.extra_args)))

subprocess.run(command, check=True, cwd=workdir, env=env)

candidates = sorted(workdir.glob("**/consensi.fa.classified"))
if not candidates:
    raise FileNotFoundError(f"RepeatModeler completed but no consensi.fa.classified was found in {workdir}")

best = max(candidates, key=lambda path: path.stat().st_size)
shutil.copy2(best, output_library)

families = 0
bases = 0
current = []
with output_library.open() as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if line.startswith(">"):
            if current:
                bases += sum(len(chunk) for chunk in current)
            families += 1
            current = []
        elif line:
            current.append(line)
if current:
    bases += sum(len(chunk) for chunk in current)

output_summary.parent.mkdir(parents=True, exist_ok=True)
output_summary.write_text(
    "sample\tlibrary_path\tfamilies\tlibrary_bases\trepeatmodeler_workdir\n"
    f"{snakemake.wildcards.accession}\t{output_library}\t{families}\t{bases}\t{workdir}\n"
)
