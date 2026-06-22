from .evo2_model import load_evo2_model
from .perbase import per_token_logprob, extract_center_window
from .export import read_fasta, write_perbase_logprobs_tsv


def export_perbase_logprobs(
    fasta_path,
    output_path="perbase_logprobs.tsv",
    model_name="evo2_7b",
    device="cuda:0",
    center=4096,
    half_window=320,
):
    model, device = load_evo2_model(model_name=model_name, device=device)

    records = read_fasta(fasta_path)
    all_rows = []

    for tag, seq in records.items():
        base_tag, allele = tag.rsplit("__", 1)

        per_token_results = per_token_logprob(
            model=model,
            seq=seq,
            device=device,
        )

        window_rows = extract_center_window(
            per_token_results,
            center=center,
            half_window=half_window,
        )

        for row in window_rows:
            row["record_tag"] = base_tag
            row["allele"] = allele

        all_rows.extend(window_rows)
        print("done", tag)

    write_perbase_logprobs_tsv(all_rows, output_path)

    print(f"Wrote {output_path}")

    return output_path