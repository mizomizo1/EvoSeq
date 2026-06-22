import re


def contains_invalid_bases(seq):
    return bool(re.search(r"[^ACGTN]", seq))


def validate_manifest_fasta_relationship(manifest_df, evo_df, strict=True):
    manifest_ids = set(manifest_df["record_id"])
    evo_ids = set(evo_df["id"])

    manifest_only = manifest_ids - evo_ids
    fasta_only = evo_ids - manifest_ids

    if fasta_only and strict:
        raise ValueError(f"Some FASTA IDs are not in manifest: {list(sorted(fasta_only))[:5]}")

    if manifest_only and strict:
        raise ValueError(f"Some manifest record_id values are not in FASTA: {list(sorted(manifest_only))[:5]}")

    if fasta_only:
        print(f"Warning: {len(fasta_only)} FASTA IDs are not in manifest.")

    if manifest_only:
        print(f"Warning: {len(manifest_only)} manifest record_id values are not in FASTA.")

    return True


def validate_sequences(evo_input_df):
    df = evo_input_df.copy()

    if df["ref_seq"].isna().any():
        raise ValueError("Missing reference sequence.")

    if df["mut_seq"].isna().any():
        raise ValueError("Missing mutant sequence.")

    df["ref_has_invalid_base"] = df["ref_seq"].apply(contains_invalid_bases)
    df["mut_has_invalid_base"] = df["mut_seq"].apply(contains_invalid_bases)

    if df["ref_has_invalid_base"].any():
        raise ValueError("Invalid base found in reference sequences.")

    if df["mut_has_invalid_base"].any():
        raise ValueError("Invalid base found in mutant sequences.")

    return df
