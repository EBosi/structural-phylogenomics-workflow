import csv
import runpy
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path("/home/bosi/kmer_phylo_workflow")
SCRIPTS_DIR = REPO_ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


class Dummy:
    pass


def run_script(script_name, snakemake_obj):
    runpy.run_path(str(SCRIPTS_DIR / script_name), init_globals={"snakemake": snakemake_obj})


def test_resolve_accessions_builds_metadata_tables(tmp_path):
    accession_file = tmp_path / "accessions.txt"
    accession_file.write_text("GCA_000001.1\nGCF_000002.1\n")
    local_genomes = tmp_path / "local_genomes.tsv"
    local_genomes.write_text("accession\tlocal_path\n")

    header = (
        "# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\t"
        "taxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\t"
        "assembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\t"
        "gbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\t"
        "relation_to_type_material\n"
    )
    genbank = tmp_path / "assembly_summary_genbank.txt"
    genbank.write_text(
        header
        + "GCA_000001.1\tPRJ1\tSAMN1\tna\tna\t11\t111\tSpecies one\tna\tna\tlatest\tScaffold\tMajor\tFull\t2024-01-01\tASM1\tLabA\tna\tna\thttps://example.org/GCA_000001.1_ASM1\tna\tna\n"
    )
    refseq = tmp_path / "assembly_summary_refseq.txt"
    refseq.write_text(
        header
        + "GCF_000002.1\tPRJ2\tSAMN2\tna\trepresentative genome\t22\t222\tSpecies two\tna\tna\tlatest\tChromosome\tMajor\tFull\t2024-02-01\tASM2\tLabB\tna\tna\thttps://example.org/GCF_000002.1_ASM2\tna\tna\n"
    )

    snk = Dummy()
    snk.input = type(
        "I",
        (),
        {
            "accession_file": str(accession_file),
            "local_genomes": str(local_genomes),
            "genbank": str(genbank),
            "refseq": str(refseq),
        },
    )()
    snk.output = type(
        "O",
        (),
        {
            "assemblies": str(tmp_path / "assemblies.tsv"),
            "organisms": str(tmp_path / "organisms.tsv"),
            "manifest": str(tmp_path / "manifest.tsv"),
        },
    )()
    snk.params = type("P", (), {"accessions": ["GCA_000001.1", "GCF_000002.1"]})()

    run_script("resolve_accessions.py", snk)

    assemblies = (tmp_path / "assemblies.tsv").read_text().splitlines()
    manifest = (tmp_path / "manifest.tsv").read_text().splitlines()
    assert len(assemblies) == 3
    assert "GCA_000001.1" in assemblies[1]
    assert "GCF_000002.1" in assemblies[2]
    assert "GCA_000001.1_ASM1_genomic.fna.gz" in manifest[1]


def test_resolve_accessions_merges_local_genome_records(tmp_path):
    local_fasta = tmp_path / "local_sample_001.fna.gz"
    local_fasta.write_text("placeholder")
    local_genomes = tmp_path / "local_genomes.tsv"
    local_genomes.write_text(
        "accession\torganism_name\tassembly_name\tassembly_level\tsource_db\tlocal_path\n"
        f"local_sample_001\tLocal species A\tLocalAsm\tScaffold\tlocal\t{local_fasta}\n"
    )

    snk = Dummy()
    snk.input = type("I", (), {"local_genomes": str(local_genomes)})()
    snk.output = type(
        "O",
        (),
        {
            "assemblies": str(tmp_path / "assemblies.tsv"),
            "organisms": str(tmp_path / "organisms.tsv"),
            "manifest": str(tmp_path / "manifest.tsv"),
        },
    )()
    snk.params = type("P", (), {"accessions": [], "eutils_base": "", "request_timeout": 60})()

    run_script("resolve_accessions.py", snk)

    assemblies = (tmp_path / "assemblies.tsv").read_text()
    manifest = (tmp_path / "manifest.tsv").read_text()
    assert "local_sample_001" in assemblies
    assert "\tlocal\t" in assemblies
    assert str(local_fasta) in manifest


def test_small_k_pipeline_scripts_produce_distance_and_tree(tmp_path):
    fasta_a = tmp_path / "A.fa"
    fasta_b = tmp_path / "B.fa"
    fasta_a.write_text(">s1\nACGTACGTACGTACGT\n")
    fasta_b.write_text(">s1\nACGTACGTAAAATTTT\n")

    for accession, fasta in (("A", fasta_a), ("B", fasta_b)):
        snk = Dummy()
        snk.input = [str(fasta)]
        snk.output = [str(tmp_path / f"{accession}.spectrum.tsv")]
        snk.params = type(
            "P",
            (),
            {
                "accession": accession,
                "dataset": "unmasked",
                "k": "3",
                "canonical": True,
                "normalization": "frequency",
            },
        )()
        run_script("compute_kmer_spectrum.py", snk)

    snk = Dummy()
    snk.input = [str(tmp_path / "A.spectrum.tsv"), str(tmp_path / "B.spectrum.tsv")]
    snk.output = [str(tmp_path / "matrix.tsv")]
    snk.params = type("P", (), {"accessions": ["A", "B"], "value_column": "frequency"})()
    run_script("build_kmer_matrix.py", snk)

    snk = Dummy()
    snk.input = [str(tmp_path / "matrix.tsv")]
    snk.output = [str(tmp_path / "distance.tsv")]
    snk.params = type("P", (), {"metric": "cosine"})()
    run_script("compute_distance_matrix.py", snk)

    snk = Dummy()
    snk.input = [str(tmp_path / "distance.tsv")]
    snk.output = [str(tmp_path / "tree.nwk")]
    snk.params = type("P", (), {"method": "nj"})()
    run_script("infer_tree.py", snk)

    distance_rows = (tmp_path / "distance.tsv").read_text().splitlines()
    tree_text = (tmp_path / "tree.nwk").read_text().strip()
    assert distance_rows[0] == "accession\tA\tB"
    assert tree_text.endswith(";")
    assert "A" in tree_text and "B" in tree_text


def test_partition_sketch_and_resampling_scripts_smoke(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">c1\nACGTACGTACGTACGT\n>c2\nACGTACGTAAAAACGT\n")

    snk = Dummy()
    snk.input = [str(fasta)]
    snk.output = type(
        "O",
        (),
        {"matrix": str(tmp_path / "windows.npy"), "meta": str(tmp_path / "windows.meta.tsv")},
    )()
    snk.params = type(
        "P",
        (),
        {
            "accession": "ACC",
            "dataset": "masked",
            "k": "3",
            "unit_type": "window",
            "window_size": 8,
            "step_size": 8,
            "min_window_length": 4,
            "canonical": True,
        },
    )()
    run_script("compute_partitioned_kmer_spectra.py", snk)
    window_matrix = np.load(tmp_path / "windows.npy", mmap_mode="r")
    assert window_matrix.shape[0] > 0

    for accession, seq in (("A", ">s\nACGTACGTACGT\n"), ("B", ">s\nACGTACGTAAAA\n")):
        path = tmp_path / f"{accession}.fa"
        path.write_text(seq)
        snk = Dummy()
        snk.input = [str(path)]
        snk.output = [str(tmp_path / f"{accession}.sig.tsv")]
        snk.params = type(
            "P",
            (),
            {"accession": accession, "dataset": "masked", "k": "5", "num_hashes": 10, "canonical": True},
        )()
        run_script("compute_minhash_signature.py", snk)

    snk = Dummy()
    snk.input = [str(tmp_path / "A.sig.tsv"), str(tmp_path / "B.sig.tsv")]
    snk.output = [str(tmp_path / "sketch.tsv")]
    snk.params = type("P", (), {"accessions": ["A", "B"]})()
    run_script("compute_sketch_distance_matrix.py", snk)
    assert "accession\tA\tB" in (tmp_path / "sketch.tsv").read_text().splitlines()[0]

    (tmp_path / "ref.nwk").write_text("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);\n")
    specs = {
        "A": [[0.8, 0.2], [0.7, 0.3]],
        "B": [[0.75, 0.25], [0.72, 0.28]],
        "C": [[0.9, 0.1], [0.85, 0.15]],
        "D": [[0.88, 0.12], [0.83, 0.17]],
    }
    matrix_paths = []
    meta_paths = []
    for accession, matrix in specs.items():
        matrix_path = tmp_path / f"{accession}.units.npy"
        meta_path = tmp_path / f"{accession}.units.meta.tsv"
        np.save(matrix_path, np.array(matrix, dtype=np.float32))
        with meta_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["accession", "dataset", "k", "unit_type", "n_units", "n_features", "dtype"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "accession": accession,
                    "dataset": "masked",
                    "k": 3,
                    "unit_type": "window",
                    "n_units": 2,
                    "n_features": 2,
                    "dtype": "float32",
                }
            )
        matrix_paths.append(str(matrix_path))
        meta_paths.append(str(meta_path))

    snk = Dummy()
    snk.input = type("I", (), {"tree": str(tmp_path / "ref.nwk"), "matrices": matrix_paths, "metadata": meta_paths})()
    snk.output = type(
        "O",
        (),
        {"support": str(tmp_path / "support.tsv"), "summary": str(tmp_path / "summary.tsv")},
    )()
    snk.params = type(
        "P",
        (),
        {"mode": "bootstrap", "metric": "cosine", "method": "nj", "replicates": 5, "fraction": 1.0, "seed": 7},
    )()
    snk.wildcards = type("W", (), {"dataset": "masked", "k": "3"})()
    run_script("resample_tree_support.py", snk)

    summary_lines = (tmp_path / "summary.tsv").read_text().splitlines()
    assert summary_lines[0].startswith("mode\tdataset\tk")
    assert "bootstrap\tmasked\t3" in summary_lines[1]


def test_sketch_distance_keeps_empty_signatures_and_skips_failed_read_signatures(tmp_path):
    empty_sig = tmp_path / "A.sig.tsv"
    empty_sig.write_text("accession\tdataset\tk\thash_rank\thash_value\n")
    populated_sig = tmp_path / "B.sig.tsv"
    populated_sig.write_text("accession\tdataset\tk\thash_rank\thash_value\nB\tmasked\t5\t1\t42\n")

    snk = Dummy()
    snk.input = [str(empty_sig), str(populated_sig)]
    snk.output = [str(tmp_path / "assembly_sketch.tsv")]
    snk.params = type("P", (), {"accessions": ["A", "B"]})()
    run_script("compute_sketch_distance_matrix.py", snk)

    rows = (tmp_path / "assembly_sketch.tsv").read_text().splitlines()
    assert rows[0] == "accession\tA\tB"
    assert rows[1] == "A\t0.000000000000\t1.000000000000"

    excluded_sig = tmp_path / "C.sig.tsv"
    excluded_sig.write_text(
        "accession\tdataset\tk\thash_rank\thash_value\tsignature_status\n"
        "C\treads_filtered\t17\t0\t\texcluded_by_qc\n"
    )
    retained_sig = tmp_path / "D.sig.tsv"
    retained_sig.write_text(
        "accession\tdataset\tk\thash_rank\thash_value\tsignature_status\n"
        "D\treads_filtered\t17\t1\t123\tok\n"
    )

    snk = Dummy()
    snk.input = [str(excluded_sig), str(retained_sig)]
    snk.output = [str(tmp_path / "reads_sketch.tsv")]
    snk.params = type(
        "P",
        (),
        {"samples": ["C", "D"], "exclude_signature_statuses": ["excluded_by_qc"]},
    )()
    run_script("compute_sketch_distance_matrix.py", snk)

    assert (tmp_path / "reads_sketch.tsv").read_text().splitlines()[0] == "accession\tD"


def test_infer_read_abundance_band_uses_dataset_k(tmp_path):
    histogram = tmp_path / "S.hist.tsv"
    histogram.write_text(
        "sample_id\tk\tsampled_reads\ttotal_kmers\tdistinct_kmers\tcount_bin\tn_distinct_kmers_at_count\n"
        "S\t15\t100\t1000\t20\t1\t10\n"
        "S\t15\t100\t1000\t20\t2\t8\n"
        "S\t15\t100\t1000\t20\t8\t2\n"
        "S\t17\t100\t1000\t63\t1\t10\n"
        "S\t17\t100\t1000\t63\t2\t2\n"
        "S\t17\t100\t1000\t63\t5\t30\n"
        "S\t17\t100\t1000\t63\t6\t20\n"
        "S\t17\t100\t1000\t63\t20\t1\n"
    )
    model = tmp_path / "S.model.tsv"
    model.write_text(
        "sample_id\tmax_supported_k\trecommended_action\n"
        "S\t15\tinclude\n"
    )
    dataset_k = tmp_path / "dataset_k.tsv"
    dataset_k.write_text(
        "n_samples\tn_included_samples\tn_excluded_samples\trecommended_k\trequired_support\t"
        "supporting_samples\tdataset_confidence\tselection_status\tselection_method\n"
        "1\t1\t0\t17\t1\t1\thigh\tselected\tquality_filtered_histogram_support\n"
    )

    snk = Dummy()
    snk.input = type(
        "I",
        (),
        {"histogram": str(histogram), "model": str(model), "dataset_k": str(dataset_k)},
    )()
    snk.output = [str(tmp_path / "S.band.tsv")]
    run_script("infer_read_abundance_band.py", snk)

    band_rows = list(csv.DictReader((tmp_path / "S.band.tsv").open(), delimiter="\t"))
    assert band_rows[0]["selected_k"] == "17"
    assert band_rows[0]["sample_supports_dataset_k"] == "0"
    assert band_rows[0]["band_confidence"] == "unsupported_dataset_k"


def test_filtered_read_signature_writes_status_for_excluded_sample(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tspecies\tread1\tread2\testimated_genome_size_bp\tplatform\tnotes\n"
        f"S\tSpecies\t{tmp_path / 'missing.fastq'}\t\t1000\tillumina\t\n"
    )
    band = tmp_path / "S.band.tsv"
    band.write_text(
        "sample_id\tselected_k\tlow_count\thigh_count\tsignal_peak_abundance\t"
        "signal_peak_count_of_counts\tsingleton_fraction\tretained_fraction\thigh_copy_fraction\t"
        "sample_supports_dataset_k\trecommended_action\tmodel_warning\tband_confidence\n"
        "S\t17\t2\t10\t0\t0\t1.0\t0.0\t0.0\t0\texclude\tno_signal_peak\texcluded\n"
    )
    dataset_k = tmp_path / "dataset_k.tsv"
    dataset_k.write_text(
        "n_samples\tn_included_samples\tn_excluded_samples\trecommended_k\trequired_support\t"
        "supporting_samples\tdataset_confidence\tselection_status\tselection_method\n"
        "1\t0\t1\t17\t0\t0\tlow\tno_included_samples\tquality_filtered_histogram_support\n"
    )

    snk = Dummy()
    snk.input = type(
        "I",
        (),
        {"manifest": str(manifest), "band": str(band), "dataset_k": str(dataset_k)},
    )()
    snk.output = [str(tmp_path / "S.signature.tsv")]
    snk.params = type("P", (), {"num_hashes": 10, "max_reads_per_file": 100, "exclude_low_confidence_bands": True})()
    snk.wildcards = type("W", (), {"sample": "S"})()
    run_script("compute_filtered_read_signature.py", snk)

    rows = list(csv.DictReader((tmp_path / "S.signature.tsv").open(), delimiter="\t"))
    assert rows[0]["signature_status"] == "excluded_by_qc"
    assert rows[0]["hash_value"] == ""


def test_repeatmasker_parser_preserves_repeat_classes(tmp_path):
    from repeat_utils import parse_repeatmasker_out_records, records_to_interval_map

    out_path = tmp_path / "repeatmasker_parser_test.out"
    out_path.write_text(
        "   SW  perc perc perc  query       position in query     matching repeat           position in repeat\n"
        "score  div. del. ins.  sequence    begin end (left)      repeat          class/family begin end (left) ID\n"
        "  500  1.2  0.0  0.0  contig1        10   40 (100) +    TE1             DNA/TcMar    1   31 (0)  1\n"
        "  300  5.0  0.0  0.0  contig1        35   60 (80)  C    Simple_repeat  Simple_repeat 1 26 (0)  2\n"
    )

    records = parse_repeatmasker_out_records(out_path, source="repeatmasker_known")
    assert records[0]["repeat_class"] == "DNA"
    assert records[0]["repeat_family"] == "TcMar"
    assert records[1]["repeat_class"] == "Simple_repeat"
    assert records_to_interval_map(records)["contig1"] == [(10, 60)]


def test_repeat_annotation_summary_reports_sources_and_classes(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">contig1\n" + "A" * 100 + "\n")
    intervals = tmp_path / "intervals.txt"
    intervals.write_text(">contig1\n1 - 30\n50 - 60\n")
    details = tmp_path / "details.tsv"
    details.write_text(
        "seq_id\tstart\tend\tsource\trepeat_name\trepeat_class\trepeat_family\tstrand\tscore\tdivergence\traw_class_family\n"
        "contig1\t1\t20\tdustmasker\tdustmasker\tLow_complexity\t\t\t\t\tLow_complexity\n"
        "contig1\t10\t30\trepeatmasker_known\tTE1\tDNA\tTcMar\t+\t500\t1.2\tDNA/TcMar\n"
        "contig1\t50\t60\trepeatmasker_denovo\tdenovo1\tLINE\tL1\t+\t400\t2.0\tLINE/L1\n"
    )

    snk = Dummy()
    snk.input = type("I", (), {"fasta": str(fasta), "intervals": str(intervals), "details": str(details)})()
    snk.output = type(
        "O",
        (),
        {"summary": str(tmp_path / "summary.tsv"), "classes": str(tmp_path / "classes.tsv")},
    )()
    snk.params = type("P", (), {"sample": "ACC"})()
    run_script("summarize_repeat_annotation.py", snk)

    summary = list(csv.DictReader((tmp_path / "summary.tsv").open(), delimiter="\t"))[0]
    classes = list(csv.DictReader((tmp_path / "classes.tsv").open(), delimiter="\t"))
    assert summary["masked_bases"] == "41"
    assert summary["dustmasker_bases"] == "20"
    assert summary["repeatmasker_known_bases"] == "21"
    assert summary["repeatmasker_denovo_bases"] == "11"
    assert {row["repeat_class"] for row in classes} == {"Low_complexity", "DNA", "LINE"}


def test_filter_organelles_removes_only_confident_contigs(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">nuc\nACGTACGT\n>mito\nTTTTCCCC\n>amb\nGGGGAAAA\n")

    calls = tmp_path / "calls.tsv"
    calls.write_text(
        "\t".join(
            [
                "query_id",
                "query_length",
                "subject_id",
                "pident",
                "aligned_length",
                "query_coverage_percent",
                "evalue",
                "bitscore",
                "classification",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "nuc\t8\t\t\t0\t0\t\t\tnuclear_like",
                "mito\t8\tref1\t99.0\t8\t100.0\t0.0\t42.0\torganelle_confident",
                "amb\t8\tref2\t82.0\t4\t50.0\t1e-20\t20.0\torganelle_ambiguous",
            ]
        )
        + "\n"
    )

    snk = Dummy()
    snk.input = type("I", (), {"fasta": str(fasta), "calls": str(calls)})()
    snk.output = type(
        "O",
        (),
        {"filtered": str(tmp_path / "filtered.fa"), "summary": str(tmp_path / "summary.tsv")},
    )()
    snk.params = type("P", (), {"accession": "ACC", "remove_classes": ["organelle_confident"]})()
    run_script("filter_organelles.py", snk)

    filtered_text = (tmp_path / "filtered.fa").read_text()
    summary_rows = (tmp_path / "summary.tsv").read_text().splitlines()

    assert ">nuc" in filtered_text
    assert ">amb" in filtered_text
    assert ">mito" not in filtered_text
    assert summary_rows[1].startswith("ACC\t2\t1\t16\t8")
    assert summary_rows[1].endswith("\t1\t1")


def test_mask_fasta_from_intervals_respects_hard_and_soft_masking(tmp_path):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">seq1\nACGTACGT\n")
    intervals = tmp_path / "intervals.txt"
    intervals.write_text(">seq1\n2 - 4\n7 - 8\n")

    snk = Dummy()
    snk.input = type("I", (), {"fasta": str(fasta), "intervals": str(intervals)})()
    snk.output = [str(tmp_path / "hard.fa")]
    snk.params = type("P", (), {"hard_masking": True})()
    run_script("mask_fasta_from_intervals.py", snk)
    assert (tmp_path / "hard.fa").read_text().splitlines()[1] == "ANNNACNN"

    snk = Dummy()
    snk.input = type("I", (), {"fasta": str(fasta), "intervals": str(intervals)})()
    snk.output = [str(tmp_path / "soft.fa")]
    snk.params = type("P", (), {"hard_masking": False})()
    run_script("mask_fasta_from_intervals.py", snk)
    assert (tmp_path / "soft.fa").read_text().splitlines()[1] == "AcgtACgt"


def test_build_pre_kmer_summary_tracks_post_organelle_metrics(tmp_path):
    assemblies = tmp_path / "assemblies.tsv"
    assemblies.write_text(
        "accession\tresolved_accession\torganism_name\tsource_db\tassembly_level\n"
        "ACC1\tACC1\tSpecies one\tgenbank\tScaffold\n"
    )
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tn_sequences\ttotal_bases\tn50\tgc_percent\n"
        "ACC1\t10\t1000\t200\t40.0\n"
    )
    preprocessing = tmp_path / "pre.tsv"
    preprocessing.write_text(
        "sample\traw_sequences\tprocessed_sequences\tretained_sequences\tremoved_sequences\t"
        "retained_sequence_fraction\traw_bases\tprocessed_bases\tretained_bases\tremoved_bases\tretained_base_fraction\n"
        "ACC1\t10\t8\t8\t2\t0.8\t1000\t900\t900\t100\t0.9\n"
    )
    organelle = tmp_path / "org.tsv"
    organelle.write_text(
        "sample\tkept_sequences\tremoved_sequences\tkept_bases\tremoved_bases\tnuclear_like_sequences\torganelle_ambiguous_sequences\torganelle_confident_sequences\n"
        "ACC1\t7\t1\t850\t50\t7\t1\t1\n"
    )
    repeats = tmp_path / "rep.tsv"
    repeats.write_text(
        "sample\tmasked_intervals\tmasked_bases\tmasked_fraction_percent\n"
        "ACC1\t12\t100\t11.7647\n"
    )

    snk = Dummy()
    snk.input = type(
        "I",
        (),
        {
            "assemblies": str(assemblies),
            "qc": str(qc),
            "preprocessing": str(preprocessing),
            "organelle": str(organelle),
            "repeats": str(repeats),
        },
    )()
    snk.output = [str(tmp_path / "pre_kmer_summary.tsv")]
    run_script("build_pre_kmer_summary.py", snk)

    rows = (tmp_path / "pre_kmer_summary.tsv").read_text().splitlines()
    assert rows[0].startswith("accession\tresolved_accession")
    assert "post_preprocess_bases" in rows[0]
    assert "post_organelle_bases" in rows[0]
    assert "\t900\t0.8\t0.9\t7\t850\t0.7\t0.85\t1\t50\t1\t1\t12\t100\t11.7647" in rows[1]
