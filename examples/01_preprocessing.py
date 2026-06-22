# examples/basic_preprocessing.py

from evoseq.preprocess import prepare_evo2_input

prepare_evo2_input(
    manifest_path="data/manifest.tsv",
    reference_fasta_path="data/evo2_reference.fasta",
    mutant_fasta_path="data/evo2_mutant.fasta",
    out_dir="evo2_input",
)