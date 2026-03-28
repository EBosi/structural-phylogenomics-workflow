import csv
import runpy
import sys
from pathlib import Path


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
        {"accession_file": str(accession_file), "genbank": str(genbank), "refseq": str(refseq)},
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
    snk.output = [str(tmp_path / "windows.tsv")]
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
    assert "unit_id" in (tmp_path / "windows.tsv").read_text().splitlines()[0]

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
    header = ["accession", "dataset", "k", "unit_type", "unit_id", "kmer", "count", "frequency"]
    specs = {
        "A": {"u1": {"AAA": 0.8, "AAC": 0.2}, "u2": {"AAA": 0.7, "AAC": 0.3}},
        "B": {"u1": {"AAA": 0.75, "AAC": 0.25}, "u2": {"AAA": 0.72, "AAC": 0.28}},
        "C": {"u1": {"CCC": 0.9, "CCG": 0.1}, "u2": {"CCC": 0.85, "CCG": 0.15}},
        "D": {"u1": {"CCC": 0.88, "CCG": 0.12}, "u2": {"CCC": 0.83, "CCG": 0.17}},
    }
    spectrum_paths = []
    for accession, units in specs.items():
        path = tmp_path / f"{accession}.units.tsv"
        spectrum_paths.append(str(path))
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for unit_id, kmers in units.items():
                for kmer, frequency in kmers.items():
                    writer.writerow(
                        {
                            "accession": accession,
                            "dataset": "masked",
                            "k": 3,
                            "unit_type": "window",
                            "unit_id": unit_id,
                            "kmer": kmer,
                            "count": 1,
                            "frequency": frequency,
                        }
                    )

    snk = Dummy()
    snk.input = type("I", (), {"tree": str(tmp_path / "ref.nwk"), "spectra": spectrum_paths})()
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
