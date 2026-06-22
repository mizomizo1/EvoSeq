from evoseq.scoring import export_perbase_logprobs

export_perbase_logprobs(
    fasta_path="representative_perbase.fasta",
    output_path="perbase_logprobs.tsv",
    model_name="evo2_7b",
    center=4096,
    half_window=320,
)
