from collections import defaultdict
from pathlib import Path

from io_utils import read_fasta


def add_interval(intervals_by_seq, seq_id, start, end):
    if start > end:
        start, end = end, start
    if start < 1:
        start = 1
    if end < 1:
        return
    intervals_by_seq[seq_id].append((start, end))


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_interval_text(path):
    intervals_by_seq = defaultdict(list)
    current_seq = None

    with Path(path).open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq = line[1:].split()[0]
                continue
            if current_seq is None or "-" not in line:
                continue
            start_str, end_str = [part.strip() for part in line.split("-", maxsplit=1)]
            add_interval(intervals_by_seq, current_seq, int(start_str), int(end_str))

    return {seq_id: merge_intervals(intervals) for seq_id, intervals in intervals_by_seq.items()}


def parse_repeatmasker_out(path):
    intervals_by_seq = defaultdict(list)

    with Path(path).open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("SW") or line.startswith("score") or line.startswith("There were no repetitive sequences"):
                continue
            parts = line.split()
            if len(parts) < 7 or not parts[0].isdigit():
                continue
            seq_id = parts[4]
            start = int(parts[5])
            end = int(parts[6])
            add_interval(intervals_by_seq, seq_id, start, end)

    return {seq_id: merge_intervals(intervals) for seq_id, intervals in intervals_by_seq.items()}


def combine_interval_sets(interval_sets):
    combined = defaultdict(list)
    for interval_map in interval_sets:
        for seq_id, intervals in interval_map.items():
            combined[seq_id].extend(intervals)
    return {seq_id: merge_intervals(intervals) for seq_id, intervals in combined.items()}


def write_interval_text(intervals_by_seq, path):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for seq_id in sorted(intervals_by_seq):
            merged = merge_intervals(intervals_by_seq[seq_id])
            if not merged:
                continue
            handle.write(f">{seq_id}\n")
            for start, end in merged:
                handle.write(f"{start} - {end}\n")


def mask_sequence(sequence, intervals, hard_masking):
    if not intervals:
        return sequence
    seq_chars = list(sequence)
    for start, end in intervals:
        start_idx = max(0, start - 1)
        end_idx = min(len(seq_chars), end)
        for idx in range(start_idx, end_idx):
            seq_chars[idx] = "N" if hard_masking else seq_chars[idx].lower()
    return "".join(seq_chars)


def write_masked_fasta(input_fasta, output_fasta, intervals_by_seq, hard_masking):
    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for header, sequence in read_fasta(input_fasta):
            seq_id = header.split()[0]
            masked = mask_sequence(sequence, intervals_by_seq.get(seq_id, []), hard_masking)
            handle.write(f">{header}\n")
            for start in range(0, len(masked), 80):
                handle.write(masked[start : start + 80] + "\n")
