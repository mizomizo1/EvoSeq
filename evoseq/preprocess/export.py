import pandas as pd

from .fasta import wrap_sequence


def sanitize_header_value(value):
    if pd.isna(value):
        return "NA"

    return str(value).replace(" ", "_").replace("\n", "_").replace("\t", "_")


def get_row_value(row, name, default="NA"):
    return getattr(row, name) if hasattr(row, name) else default


def make_fasta_header(row, allele):
    seq_len = row.ref_len if allele == "ref" else row.mut_len

    fields = [
        row.id,
        f"allele={allele}",
        f"gene={sanitize_header_value(get_row_value(row, 'gene'))}",
        f"variant={sanitize_header_value(get_row_value(row, 'variant'))}",
        f"hgvs={sanitize_header_value(get_row_value(row, 'hgvs'))}",
        f"ann={sanitize_header_value(get_row_value(row, 'annotation'))}",
        f"type={sanitize_header_value(get_row_value(row, 'variant_type'))}",
        f"len={seq_len}",
    ]

    return "|".join(fields)


def write_fasta_from_df(table, path, allele):
    seq_col = "ref_seq" if allele == "ref" else "mut_seq"

    with open(path, "w") as f:
        for row in table.itertuples(index=False):
            header = make_fasta_header(row, allele)
            seq = getattr(row, seq_col)

            f.write(f">{header}\n")
            f.write(wrap_sequence(seq) + "\n")


def export_evo2_input(evo_input_df, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs_path = out_dir / "evo2_pairs.tsv"
    pair_path = out_dir / "evo2_pair.tsv"
    ref_path = out_dir / "evo2_reference.fa"
    mut_path = out_dir / "evo2_mutant.fa"
    all_path = out_dir / "evo2_all.fa"

    evo_input_df.to_csv(pairs_path, sep="\t", index=False)
    evo_input_df.to_csv(pair_path, sep="\t", index=False)

    write_fasta_from_df(evo_input_df, ref_path, "ref")
    write_fasta_from_df(evo_input_df, mut_path, "mut")

    with open(all_path, "w") as fout:
        for row in evo_input_df.itertuples(index=False):
            ref_header = make_fasta_header(row, "ref")
            mut_header = make_fasta_header(row, "mut")

            fout.write(f">{ref_header}\n")
            fout.write(wrap_sequence(row.ref_seq) + "\n")
            fout.write(f">{mut_header}\n")
            fout.write(wrap_sequence(row.mut_seq) + "\n")

    return {
        "pairs": pairs_path,
        "pair_tsv": pair_path,
        "reference": ref_path,
        "reference_fasta": ref_path,
        "mutant": mut_path,
        "mutant_fasta": mut_path,
        "all": all_path,
        "all_fasta": all_path,
    }
