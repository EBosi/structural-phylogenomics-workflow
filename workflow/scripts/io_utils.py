import gzip


def open_maybe_gzip(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    name = None
    seq_chunks = []
    with open_maybe_gzip(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if name is not None:
            yield name, "".join(seq_chunks)


def read_fastq(path):
    record_number = 0
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            record_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record {record_number} in {path}")
            header = header.rstrip()
            sequence = sequence.rstrip()
            plus = plus.rstrip()
            quality = quality.rstrip()
            if not header.startswith("@"):
                raise ValueError(f"FASTQ record {record_number} in {path} has an invalid header")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ record {record_number} in {path} has an invalid plus line")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"FASTQ record {record_number} in {path} has sequence/quality length mismatch"
                )
            yield header[1:], sequence, quality
