#!/usr/bin/env python3
"""
Compute a diagnostic metric panel on existing retained sketch hash sets.

Important:
- This reads existing signature files and treats retained hashes as sets.
- Metrics are computed on retained sketch hashes, not exact full k-mer sets.
- containment_a_in_b here means: fraction of retained hashes from sample A also present in sample B.
- This is useful as a diagnostic, but it is not rigorous full-set containment.
- If all signatures have the same number of retained hashes, many metrics become monotonic functions of the same intersection count.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from glob import glob


ACCESSION_RE = re.compile(r"(protura_purged|GC[AF]_\d+\.\d+)")

OUTPUT_COLUMNS = [
    "sample_a",
    "sample_b",
    "n_a",
    "n_b",
    "intersection",
    "union",
    "jaccard",
    "jaccard_distance",
    "dice",
    "dice_distance",
    "overlap_coefficient",
    "overlap_distance",
    "containment_a_in_b",
    "containment_b_in_a",
    "max_containment",
    "binary_cosine",
    "binary_cosine_distance",
    "mash_distance",
]

SIMILARITY_METRICS = [
    "jaccard",
    "dice",
    "overlap_coefficient",
    "containment_a_in_b",
    "containment_b_in_a",
    "max_containment",
    "binary_cosine",
]

DISTANCE_METRICS = [
    "jaccard_distance",
    "dice_distance",
    "overlap_distance",
    "binary_cosine_distance",
    "mash_distance",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--signatures-glob",
        action="append",
        required=True,
        help=(
            "Glob for signature files. Can be used multiple times. "
            "Quote it so Python expands it, e.g. "
            "'results/phase2_expanded_panel/**/*k15*masked*'"
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for diagnostic TSV files.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=15,
        help="k-mer size used for Mash distance conversion. Default: 15.",
    )
    parser.add_argument(
        "--hash-column",
        default="auto",
        help=(
            "Hash column name or zero-based column index. "
            "Default: auto. Auto prefers a column whose header contains 'hash'."
        ),
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Allow empty signatures instead of failing.",
    )

    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def split_line(line: str) -> list[str]:
    line = line.strip()
    if "\t" in line:
        return line.split("\t")
    if "," in line:
        return line.split(",")
    return line.split()


def parse_int_token(token: str) -> int | None:
    token = token.strip()
    if not token:
        return None
    try:
        if token.lower().startswith("0x"):
            return int(token, 16)
        return int(token)
    except ValueError:
        return None


def infer_sample_name(path: Path) -> str:
    """
    Prefer accession-like names anywhere in the path.
    Fallback to filename stem with common suffix cleanup.
    """
    for part in reversed(path.parts):
        match = ACCESSION_RE.search(part)
        if match:
            return match.group(1)

    name = path.name
    for suffix in [
        ".tsv.gz",
        ".txt.gz",
        ".csv.gz",
        ".sig.gz",
        ".signature.gz",
        ".tsv",
        ".txt",
        ".csv",
        ".sig",
        ".signature",
    ]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return Path(name).stem


def detect_hash_column(tokens: list[str], hash_column: str) -> tuple[int | None, bool]:
    """
    Return (column_index, first_line_is_header).
    """
    if hash_column != "auto":
        try:
            return int(hash_column), False
        except ValueError:
            wanted = hash_column.lower()
            for i, token in enumerate(tokens):
                if token.strip().lower() == wanted:
                    return i, True
            raise ValueError(f"Requested hash column not found in header: {hash_column}")

    lowered = [t.strip().lower() for t in tokens]

    for preferred in ("hash_value", "hash", "hashvalue"):
        if preferred in lowered:
            return lowered.index(preferred), True

    for i, token in enumerate(lowered):
        if "hash" in token and "rank" not in token:
            return i, True

    numeric_positions = [i for i, token in enumerate(tokens) if parse_int_token(token) is not None]
    if numeric_positions:
        return numeric_positions[0], False

    raise ValueError(
        "Could not auto-detect a hash column. "
        "Use --hash-column with a column name or zero-based index."
    )


def read_hashes(path: Path, hash_column: str) -> set[int]:
    hashes: set[int] = set()

    with open_text(path) as handle:
        first_tokens = None

        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            first_tokens = split_line(line)
            break

        if first_tokens is None:
            return hashes

        col_idx, first_line_is_header = detect_hash_column(first_tokens, hash_column)

        if col_idx is None:
            raise ValueError(f"Could not detect hash column for {path}")

        if not first_line_is_header:
            if col_idx >= len(first_tokens):
                raise ValueError(f"Column index {col_idx} out of range in {path}")
            value = parse_int_token(first_tokens[col_idx])
            if value is None:
                raise ValueError(f"First data row does not contain an integer hash in {path}")
            hashes.add(value)

        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            tokens = split_line(line)
            if col_idx >= len(tokens):
                continue
            value = parse_int_token(tokens[col_idx])
            if value is not None:
                hashes.add(value)

    return hashes


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0:
        return None
    return num / den


def mash_distance_from_jaccard(jaccard: float | None, k: int) -> float | None:
    if jaccard is None:
        return None
    if jaccard == 0:
        return math.inf
    return -(1.0 / k) * math.log((2.0 * jaccard) / (1.0 + jaccard))


def compute_metrics(set_a: set[int], set_b: set[int], k: int) -> dict[str, float | int | None]:
    n_a = len(set_a)
    n_b = len(set_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)

    jaccard = safe_ratio(intersection, union)
    dice = safe_ratio(2 * intersection, n_a + n_b)
    overlap = safe_ratio(intersection, min(n_a, n_b))
    containment_a_in_b = safe_ratio(intersection, n_a)
    containment_b_in_a = safe_ratio(intersection, n_b)
    max_containment = safe_ratio(intersection, max(n_a, n_b))
    binary_cosine = safe_ratio(intersection, math.sqrt(n_a * n_b))

    return {
        "n_a": n_a,
        "n_b": n_b,
        "intersection": intersection,
        "union": union,
        "jaccard": jaccard,
        "jaccard_distance": None if jaccard is None else 1.0 - jaccard,
        "dice": dice,
        "dice_distance": None if dice is None else 1.0 - dice,
        "overlap_coefficient": overlap,
        "overlap_distance": None if overlap is None else 1.0 - overlap,
        "containment_a_in_b": containment_a_in_b,
        "containment_b_in_a": containment_b_in_a,
        "max_containment": max_containment,
        "binary_cosine": binary_cosine,
        "binary_cosine_distance": None if binary_cosine is None else 1.0 - binary_cosine,
        "mash_distance": mash_distance_from_jaccard(jaccard, k),
    }


def format_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        if value == 0.0:
            return "0.0"
        return repr(value)
    return str(value)


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def expand_signature_paths(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(Path(p) for p in glob(pattern, recursive=True))
    paths = sorted(set(p for p in paths if p.is_file()))
    return paths


def choose_best(rows: list[dict], metric: str) -> dict | None:
    valid = [r for r in rows if r["sample_a"] != r["sample_b"] and r[metric] is not None]
    if not valid:
        return None

    if metric in DISTANCE_METRICS:
        return min(valid, key=lambda r: r[metric])
    return max(valid, key=lambda r: r[metric])


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    signature_paths = expand_signature_paths(args.signatures_glob)
    if not signature_paths:
        print("ERROR: no signature files matched --signatures-glob", file=sys.stderr)
        return 2

    samples: list[tuple[str, Path, set[int]]] = []
    seen_names: dict[str, Path] = {}

    for path in signature_paths:
        sample = infer_sample_name(path)
        if sample in seen_names:
            print(
                f"ERROR: duplicate inferred sample name '{sample}' for:\n"
                f"  {seen_names[sample]}\n"
                f"  {path}\n"
                "Use a narrower glob or rename/copy files into sample-specific paths.",
                file=sys.stderr,
            )
            return 2

        hashes = read_hashes(path, args.hash_column)
        if not hashes and not args.allow_empty:
            print(f"ERROR: empty signature: {path}", file=sys.stderr)
            return 2

        seen_names[sample] = path
        samples.append((sample, path, hashes))

    size_rows = [
        {
            "sample": sample,
            "signature_path": str(path),
            "n_hashes": len(hashes),
        }
        for sample, path, hashes in samples
    ]
    write_tsv(
        outdir / "signature_size_summary.tsv",
        size_rows,
        ["sample", "signature_path", "n_hashes"],
    )

    pair_rows_raw: list[dict] = []
    pair_rows_formatted: list[dict] = []

    for sample_a, _, set_a in samples:
        for sample_b, _, set_b in samples:
            metrics = compute_metrics(set_a, set_b, args.k)
            raw_row = {"sample_a": sample_a, "sample_b": sample_b}
            raw_row.update(metrics)
            pair_rows_raw.append(raw_row)

            formatted_row = {"sample_a": sample_a, "sample_b": sample_b}
            formatted_row.update({k: format_value(v) for k, v in metrics.items()})
            pair_rows_formatted.append(formatted_row)

    write_tsv(outdir / "sketch_metric_panel.tsv", pair_rows_formatted, OUTPUT_COLUMNS)

    nearest_rows = []
    metrics_to_rank = SIMILARITY_METRICS + DISTANCE_METRICS

    by_focal: dict[str, list[dict]] = defaultdict(list)
    for row in pair_rows_raw:
        by_focal[row["sample_a"]].append(row)

    for focal, rows in sorted(by_focal.items()):
        for metric in metrics_to_rank:
            best = choose_best(rows, metric)
            if best is None:
                continue
            nearest_rows.append(
                {
                    "focal_sample": focal,
                    "metric": metric,
                    "best_match": best["sample_b"],
                    "score": format_value(best[metric]),
                    "rank_rule": "min" if metric in DISTANCE_METRICS else "max",
                }
            )

    write_tsv(
        outdir / "nearest_by_metric.tsv",
        nearest_rows,
        ["focal_sample", "metric", "best_match", "score", "rank_rule"],
    )

    sizes = [len(hashes) for _, _, hashes in samples]
    notes = []
    notes.append(f"samples\t{len(samples)}")
    notes.append(f"min_hashes\t{min(sizes)}")
    notes.append(f"max_hashes\t{max(sizes)}")
    notes.append(f"all_equal_size\t{str(len(set(sizes)) == 1).lower()}")
    notes.append("")
    notes.append(
        "IMPORTANT: metrics are computed on retained sketch hashes, not exact full k-mer sets."
    )
    notes.append(
        "IMPORTANT: containment values here are sketch-set diagnostics, not rigorous full-set containment estimates."
    )

    if len(set(sizes)) == 1:
        notes.append("")
        notes.append(
            "WARNING: all signatures have the same retained hash count. "
            "In this case dice, overlap, directional containment, max_containment, "
            "and binary_cosine are mostly monotonic functions of the same intersection count. "
            "They may not add independent ranking information."
        )
    else:
        notes.append("")
        notes.append(
            "NOTE: signature sizes vary. Directional containment and overlap may differ from symmetric metrics, "
            "but interpretation is still limited by the retained-sketch design."
        )

    (outdir / "diagnostic_notes.txt").write_text("\n".join(notes) + "\n")

    print(f"Wrote: {outdir / 'signature_size_summary.tsv'}")
    print(f"Wrote: {outdir / 'sketch_metric_panel.tsv'}")
    print(f"Wrote: {outdir / 'nearest_by_metric.tsv'}")
    print(f"Wrote: {outdir / 'diagnostic_notes.txt'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
