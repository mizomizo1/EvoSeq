from pathlib import Path

from .manifest import read_manifest, aggregate_manifest
from .fasta import read_fasta_as_dict, check_fasta_pair
from .variant import compare_ref_mut, build_variant_table
from .validation import validate_manifest_fasta_relationship, validate_sequences
from .export import export_evo2_input


def prepare_evo2_input(
    manifest_path,
    reference_fasta_path,
    mutant_fasta_path,
    out_dir,
):
    manifest_path = Path(manifest_path)
    reference_fasta_path = Path(reference_fasta_path)
    mutant_fasta_path = Path(mutant_fasta_path)
    out_dir = Path(out_dir)

    manifest_df = read_manifest(manifest_path)

    ref_records = read_fasta_as_dict(reference_fasta_path)
    mut_records = read_fasta_as_dict(mutant_fasta_path)

    common_ids = check_fasta_pair(ref_records, mut_records)

    pair_df = compare_ref_mut(ref_records, mut_records, common_ids)
    variant_df = build_variant_table(pair_df)

    evo_df = variant_df[
        [
            "id",
            "source",
            "sample_from_id",
            "gene",
            "locus",
            "variant",
            "variant_type",
            "ref_len",
            "mut_len",
            "first_diff_0based",
            "first_diff_1based",
            "ref_base",
            "mut_base",
        ]
    ].copy()

    evo_df["ref_seq"] = evo_df["id"].map(ref_records)
    evo_df["mut_seq"] = evo_df["id"].map(mut_records)

    validate_manifest_fasta_relationship(manifest_df, evo_df)

    manifest_unique = aggregate_manifest(manifest_df)

    evo_input_df = evo_df.merge(
        manifest_unique,
        left_on="id",
        right_on="record_id",
        how="left",
    )

    evo_input_df = validate_sequences(evo_input_df)

    output_paths = export_evo2_input(evo_input_df, out_dir)

    return evo_input_df, output_paths