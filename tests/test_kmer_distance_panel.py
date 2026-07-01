import math
import subprocess
import sys
from pathlib import Path

import csv


def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def canonical_kmer(kmer):
    rc = reverse_complement(kmer)
    return min(kmer, rc)


def unique_kmer_set(sequences, k, canonical=True):
    kmers = set()
    for sequence in sequences:
        seq = sequence.upper()
        if len(seq) < k:
            continue
        for idx in range(len(seq) - k + 1):
            kmer = seq[idx : idx + k]
            if any(base not in "ACGT" for base in kmer):
                continue
            if canonical:
                kmer = canonical_kmer(kmer)
            kmers.add(kmer)
    return kmers


def safe_ratio(numerator, denominator):
    if denominator == 0:
        return None
    return numerator / denominator


def mash_distance_from_jaccard(jaccard, k):
    if jaccard is None:
        return None
    if jaccard == 0:
        return math.inf
    return -(1.0 / k) * math.log((2.0 * jaccard) / (1.0 + jaccard))


def compute_distance_panel(set_a, set_b, k):
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


def assert_zero_distances(panel):
    for key in (
        "jaccard_distance",
        "dice_distance",
        "overlap_distance",
        "binary_cosine_distance",
        "mash_distance",
    ):
        assert panel[key] == 0.0


def test_identical_sequences_have_perfect_similarity():
    k = 5
    seq = ["ACGTTGCATGTCAGTCA"]
    set_a = unique_kmer_set(seq, k)
    set_b = unique_kmer_set(seq, k)
    panel = compute_distance_panel(set_a, set_b, k)

    assert panel["intersection"] == panel["union"] == panel["n_a"] == panel["n_b"]
    for key in (
        "jaccard",
        "dice",
        "overlap_coefficient",
        "containment_a_in_b",
        "containment_b_in_a",
        "max_containment",
        "binary_cosine",
    ):
        assert panel[key] == 1.0
    assert_zero_distances(panel)


def test_disjoint_sequences_have_zero_similarity():
    k = 3
    set_a = unique_kmer_set(["AAAAAA"], k)
    set_b = unique_kmer_set(["CCCCCC"], k)
    panel = compute_distance_panel(set_a, set_b, k)

    assert panel["intersection"] == 0
    assert panel["union"] == panel["n_a"] + panel["n_b"]
    for key in (
        "jaccard",
        "dice",
        "overlap_coefficient",
        "containment_a_in_b",
        "containment_b_in_a",
        "max_containment",
        "binary_cosine",
    ):
        assert panel[key] == 0.0
    for key in (
        "jaccard_distance",
        "dice_distance",
        "overlap_distance",
        "binary_cosine_distance",
    ):
        assert panel[key] == 1.0
    assert math.isinf(panel["mash_distance"])


def test_subset_case_exposes_jaccard_containment_mismatch():
    k = 5
    seq_a = ["ACGTTGCATGTCAGTCA"]
    seq_b = ["ACGTTGCATGTCAGTCATTTGGGCCCAAATTT"]
    set_a = unique_kmer_set(seq_a, k)
    set_b = unique_kmer_set(seq_b, k)

    assert set_a < set_b
    panel = compute_distance_panel(set_a, set_b, k)

    assert panel["containment_a_in_b"] == 1.0
    assert panel["overlap_coefficient"] == 1.0
    assert panel["containment_b_in_a"] == panel["jaccard"] == panel["n_a"] / panel["n_b"]
    assert 0.0 < panel["dice"] < 1.0
    assert panel["jaccard"] < panel["dice"] < panel["binary_cosine"] <= panel["overlap_coefficient"]


def test_shared_core_plus_private_tails_separates_symmetric_metrics():
    k = 4
    seq_a = ["AAAACCCCGGGGAAAA"]
    seq_b = ["CCCCGGGGTATATATA"]
    set_a = unique_kmer_set(seq_a, k)
    set_b = unique_kmer_set(seq_b, k)
    panel = compute_distance_panel(set_a, set_b, k)

    assert 0 < panel["intersection"] < min(panel["n_a"], panel["n_b"])
    assert panel["jaccard"] < panel["dice"] <= panel["binary_cosine"] <= panel["overlap_coefficient"]
    assert panel["containment_a_in_b"] < 1.0
    assert panel["containment_b_in_a"] < 1.0


def test_size_imbalance_depresses_jaccard_but_not_containment_of_smaller_set():
    k = 5
    seq_a = ["ACGTTGCATGTCAGTCA"]
    seq_b_small_tail = ["ACGTTGCATGTCAGTCATTTGGGCCCAAATTT"]
    seq_b_large_tail = ["ACGTTGCATGTCAGTCATTTGGGCCCAAATTTGGGGAAAACCCCTTTTGGGGAAAACCCC"]

    set_a = unique_kmer_set(seq_a, k)
    set_b1 = unique_kmer_set(seq_b_small_tail, k)
    set_b2 = unique_kmer_set(seq_b_large_tail, k)

    assert set_a < set_b1
    assert set_a < set_b2
    panel_b1 = compute_distance_panel(set_a, set_b1, k)
    panel_b2 = compute_distance_panel(set_a, set_b2, k)

    assert panel_b1["containment_a_in_b"] == 1.0
    assert panel_b2["containment_a_in_b"] == 1.0
    assert panel_b2["n_b"] > panel_b1["n_b"]
    assert panel_b2["jaccard"] < panel_b1["jaccard"]
    assert panel_b2["dice"] < panel_b1["dice"]
    assert panel_b2["binary_cosine"] < panel_b1["binary_cosine"]


def test_fragmentation_loses_boundary_kmers_even_for_same_underlying_sequence():
    k = 5
    full = unique_kmer_set(["ACGTTGCATGTCAGTCA"], k)
    fragmented = unique_kmer_set(["ACGTTGCAT", "GTCAGTCA"], k)
    panel = compute_distance_panel(full, fragmented, k)

    assert fragmented < full
    assert panel["containment_b_in_a"] == 1.0
    assert panel["containment_a_in_b"] < 1.0
    assert panel["jaccard"] < 1.0
    assert panel["overlap_coefficient"] == 1.0


def test_repeat_copy_number_does_not_change_unique_kmer_set_metrics():
    k = 3
    set_a = unique_kmer_set(["ATGATGATG"], k)
    set_b = unique_kmer_set(["ATGATGATGATGATGATG"], k)
    panel = compute_distance_panel(set_a, set_b, k)

    assert set_a == set_b
    assert panel["n_a"] == panel["n_b"]
    for key in (
        "jaccard",
        "dice",
        "overlap_coefficient",
        "containment_a_in_b",
        "containment_b_in_a",
        "max_containment",
        "binary_cosine",
    ):
        assert panel[key] == 1.0


def test_zero_size_sets_are_handled_explicitly():
    k = 5
    empty = unique_kmer_set(["AAA"], k)
    nonempty = unique_kmer_set(["ACGTTGCATGTCAGTCA"], k)

    panel_empty_empty = compute_distance_panel(empty, empty, k)
    assert panel_empty_empty["n_a"] == 0
    assert panel_empty_empty["n_b"] == 0
    for key in (
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
    ):
        assert panel_empty_empty[key] is None

    panel_empty_nonempty = compute_distance_panel(empty, nonempty, k)
    assert panel_empty_nonempty["jaccard"] == 0.0
    assert panel_empty_nonempty["jaccard_distance"] == 1.0
    assert panel_empty_nonempty["dice"] == 0.0
    assert panel_empty_nonempty["dice_distance"] == 1.0
    assert panel_empty_nonempty["containment_a_in_b"] is None
    assert panel_empty_nonempty["containment_b_in_a"] == 0.0
    assert panel_empty_nonempty["overlap_coefficient"] is None
    assert panel_empty_nonempty["binary_cosine"] is None
    assert math.isinf(panel_empty_nonempty["mash_distance"])


def run_distance_panel_script(tmp_path, args):
    script = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "compute_kmer_set_distance_panel.py"
    output = tmp_path / "panel.tsv"
    cmd = [sys.executable, str(script), "--output", str(output), *args]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return completed, rows


def write_fasta(path, records):
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def test_script_identical_fasta_gives_zero_distance(tmp_path):
    fasta_a = tmp_path / "A.fa"
    fasta_b = tmp_path / "B.fa"
    write_fasta(fasta_a, [("a", "ACGTTGCATGTCAGTCA")])
    write_fasta(fasta_b, [("b", "ACGTTGCATGTCAGTCA")])

    _, rows = run_distance_panel_script(
        tmp_path,
        ["--k", "5", "--canonical", "true", str(fasta_a), str(fasta_b)],
    )

    assert len(rows) == 4
    pair = next(row for row in rows if row["sample_a"] == "A" and row["sample_b"] == "B")
    assert pair["jaccard_distance"] == "0.0"
    assert pair["mash_distance"] == "0.0"


def test_script_subset_case_preserves_containment_signal(tmp_path):
    fasta_a = tmp_path / "A.fa"
    fasta_b = tmp_path / "B.fa"
    write_fasta(fasta_a, [("a", "ACGTTGCATGTCAGTCA")])
    write_fasta(fasta_b, [("b", "ACGTTGCATGTCAGTCATTTGGGCCCAAATTT")])

    _, rows = run_distance_panel_script(
        tmp_path,
        ["--k", "5", "--canonical", "true", str(fasta_a), str(fasta_b)],
    )

    pair = next(row for row in rows if row["sample_a"] == "A" and row["sample_b"] == "B")
    assert pair["containment_a_in_b"] == "1.0"
    assert float(pair["jaccard"]) < 1.0


def test_script_fragmented_input_loses_boundary_kmers(tmp_path):
    fasta_full = tmp_path / "full.fa"
    fasta_frag = tmp_path / "frag.fa"
    write_fasta(fasta_full, [("full", "ACGTTGCATGTCAGTCA")])
    write_fasta(fasta_frag, [("f1", "ACGTTGCAT"), ("f2", "GTCAGTCA")])

    _, rows = run_distance_panel_script(
        tmp_path,
        ["--k", "5", "--canonical", "true", str(fasta_full), str(fasta_frag)],
    )

    pair = next(row for row in rows if row["sample_a"] == "full" and row["sample_b"] == "frag")
    assert pair["overlap_coefficient"] == "1.0"
    assert float(pair["jaccard"]) < 1.0


def test_script_skips_non_acgt_kmers(tmp_path):
    fasta_a = tmp_path / "clean.fa"
    fasta_b = tmp_path / "dirty.fa"
    write_fasta(fasta_a, [("clean", "ACGTTGCATGTCAGTCA")])
    write_fasta(fasta_b, [("dirty", "ACGNNNTGCATGTCAGTCA")])

    _, rows = run_distance_panel_script(
        tmp_path,
        ["--k", "5", "--canonical", "true", str(fasta_a), str(fasta_b)],
    )

    pair = next(row for row in rows if row["sample_a"] == "clean" and row["sample_b"] == "dirty")
    assert float(pair["intersection"]) > 0.0
    assert float(pair["jaccard"]) < 1.0


def test_script_accepts_gzipped_fasta_and_reports_expected_columns(tmp_path):
    import gzip

    fasta_a = tmp_path / "A.fa.gz"
    fasta_b = tmp_path / "B.fa.gz"
    with gzip.open(fasta_a, "wt") as handle:
        handle.write(">a\nACGTTGCATGTCAGTCA\n")
    with gzip.open(fasta_b, "wt") as handle:
        handle.write(">b\nCCCCCCCCCCCCCCCCC\n")

    _, rows = run_distance_panel_script(
        tmp_path,
        ["--k", "5", "--canonical", "false", str(fasta_a), str(fasta_b)],
    )

    assert rows
    expected = [
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
    assert list(rows[0].keys()) == expected


def test_script_fails_on_duplicate_inferred_sample_names(tmp_path):
    fasta_a = tmp_path / "dup.fa"
    fasta_b_dir = tmp_path / "nested"
    fasta_b_dir.mkdir()
    fasta_b = fasta_b_dir / "dup.fa.gz"
    write_fasta(fasta_a, [("a", "ACGTTGCATGTCAGTCA")])
    import gzip

    with gzip.open(fasta_b, "wt") as handle:
        handle.write(">b\nACGTTGCATGTCAGTCA\n")

    script = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "compute_kmer_set_distance_panel.py"
    output = tmp_path / "panel.tsv"
    cmd = [
        sys.executable,
        str(script),
        "--output",
        str(output),
        "--k",
        "5",
        "--canonical",
        "true",
        str(fasta_a),
        str(fasta_b),
    ]
    completed = subprocess.run(cmd, capture_output=True, text=True)

    assert completed.returncode != 0
    assert "Duplicate sample names" in completed.stderr
