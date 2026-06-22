from pathlib import Path

from .manifest import read_manifest, aggregate_manifest
from .fasta import read_fasta_as_dict, check_fasta_pair
from .variant import compare_ref_mut, build_variant_table
from .validation import validate_manifest_fasta_relationship, validate_sequences
from .export import export_evo2_input
from .discovery import discover_fasta_pair, discover_manifest, infer_dataset_type
from ..paths import default_output_dir, ensure_output_dir


def _progress(iterable, enabled=True, **kwargs):
    if not enabled:
        return iterable
    try:
        from tqdm.auto import tqdm

        return tqdm(iterable, **kwargs)
    except Exception:
        return iterable


def prepare_output_dir(path):
    return ensure_output_dir(path, fallback="/content/evoseq_preprocess_output")


def write_preprocessing_report(evo_input_df, output_paths, out_dir, manifest_path=None):
    report_path = Path(out_dir) / "preprocessing_report.tsv"
    report = {
        "n_variants": len(evo_input_df),
        "n_unique_ids": evo_input_df["id"].nunique(),
        "n_duplicated_ids": int(evo_input_df["id"].duplicated().sum()),
        "n_ref_invalid": int(evo_input_df["ref_has_invalid_base"].sum()),
        "n_mut_invalid": int(evo_input_df["mut_has_invalid_base"].sum()),
        "min_ref_len": int(evo_input_df["ref_len"].min()),
        "max_ref_len": int(evo_input_df["ref_len"].max()),
        "min_mut_len": int(evo_input_df["mut_len"].min()),
        "max_mut_len": int(evo_input_df["mut_len"].max()),
        "manifest_path": str(manifest_path) if manifest_path else "",
        "pairs_path": str(output_paths["pairs"]),
        "reference_fasta": str(output_paths["reference"]),
        "mutant_fasta": str(output_paths["mutant"]),
        "all_fasta": str(output_paths["all"]),
    }

    if "variant_type" in evo_input_df.columns:
        for name, count in evo_input_df["variant_type"].value_counts(dropna=False).items():
            report[f"variant_type_{name}"] = int(count)

    if "first_diff_0based" in evo_input_df.columns:
        report["median_first_diff_0based"] = float(
            evo_input_df["first_diff_0based"].dropna().median()
        )

    import pandas as pd

    pd.DataFrame([report]).to_csv(report_path, sep="\t", index=False)
    output_paths["preprocessing_report"] = report_path
    return report_path


def prepare_evo2_input(
    manifest_path=None,
    reference_fasta_path=None,
    mutant_fasta_path=None,
    out_dir=None,
    strict_manifest=True,
    progress=True,
):
    if reference_fasta_path is None or mutant_fasta_path is None:
        raise ValueError("reference_fasta_path and mutant_fasta_path are required.")

    reference_fasta_path = Path(reference_fasta_path)
    mutant_fasta_path = Path(mutant_fasta_path)
    if out_dir is None:
        out_dir = default_output_dir(
            "preprocess",
            reference_fasta_path,
            mutant_fasta_path,
            manifest_path,
        )
    out_dir = prepare_output_dir(out_dir)

    print("Preparing Evo2 input files.")
    print(f"Reference FASTA: {reference_fasta_path}")
    print(f"Mutant FASTA   : {mutant_fasta_path}")
    if manifest_path:
        manifest_path = Path(manifest_path)
        print(f"Manifest       : {manifest_path}")
        manifest_df = read_manifest(manifest_path)
    else:
        print("Manifest       : not provided; metadata will be inferred from FASTA IDs.")
        manifest_df = None

    ref_records = read_fasta_as_dict(reference_fasta_path)
    mut_records = read_fasta_as_dict(mutant_fasta_path)

    print(f"Reference records: {len(ref_records)}")
    print(f"Mutant records   : {len(mut_records)}")

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

    evo_df["ref_seq"] = evo_df["id"].map(ref_records).str.upper()
    evo_df["mut_seq"] = evo_df["id"].map(mut_records).str.upper()

    if manifest_df is not None:
        validate_manifest_fasta_relationship(
            manifest_df,
            evo_df,
            strict=strict_manifest,
        )

        manifest_unique = aggregate_manifest(manifest_df)

        evo_input_df = evo_df.merge(
            manifest_unique,
            left_on="id",
            right_on="record_id",
            how="left",
        )
    else:
        evo_input_df = evo_df

    evo_input_df = validate_sequences(evo_input_df)

    output_paths = export_evo2_input(evo_input_df, out_dir)
    report_path = write_preprocessing_report(
        evo_input_df,
        output_paths,
        out_dir,
        manifest_path=manifest_path,
    )

    print("Saved Evo2 input files:")
    for name, path in output_paths.items():
        print(f"  {name}: {path}")
    print(f"Preprocessing report: {report_path}")

    return evo_input_df, output_paths


def preprocess_files(
    reference_fasta_path,
    mutant_fasta_path,
    manifest_path="auto",
    output_dir=None,
    strict_manifest=False,
    progress=True,
):
    if manifest_path == "auto":
        manifest_path = discover_manifest(
            Path(reference_fasta_path).parent,
            manifest_path="auto",
        )

    if output_dir is None:
        output_dir = default_output_dir(
            "preprocess",
            reference_fasta_path,
            mutant_fasta_path,
            manifest_path,
        )

    print("Running EvoSeq preprocessing from explicit files.")
    print(f"Output directory: {output_dir}")

    return prepare_evo2_input(
        reference_fasta_path=reference_fasta_path,
        mutant_fasta_path=mutant_fasta_path,
        manifest_path=manifest_path,
        out_dir=output_dir,
        strict_manifest=strict_manifest,
        progress=progress,
    )


def preprocess_from_base_dir(
    base_dir,
    out_dir=None,
    manifest_path="auto",
    reference_fasta_path=None,
    mutant_fasta_path=None,
    dataset_type="auto",
    window_size=None,
    strict_manifest=False,
    progress=True,
):
    base_dir = Path(base_dir)
    if out_dir is None:
        out_dir = default_output_dir("preprocess", base_dir=base_dir)

    if dataset_type == "auto":
        dataset_type = infer_dataset_type(base_dir)

    manifest = discover_manifest(base_dir, manifest_path=manifest_path)
    reference, mutant = discover_fasta_pair(
        base_dir,
        reference_fasta_path=reference_fasta_path,
        mutant_fasta_path=mutant_fasta_path,
        dataset_type=dataset_type,
        window_size=window_size,
    )

    print(f"Base directory : {base_dir}")
    print(f"Dataset type   : {dataset_type}")

    return prepare_evo2_input(
        reference_fasta_path=reference,
        mutant_fasta_path=mutant,
        manifest_path=manifest,
        out_dir=out_dir,
        strict_manifest=strict_manifest,
        progress=progress,
    )


def preprocess_folder(
    input_dir,
    output_dir=None,
    manifest_path="auto",
    reference_fasta_path=None,
    mutant_fasta_path=None,
    dataset_type="auto",
    window_size=None,
    strict_manifest=False,
    progress=True,
):
    return preprocess_from_base_dir(
        base_dir=input_dir,
        out_dir=output_dir,
        manifest_path=manifest_path,
        reference_fasta_path=reference_fasta_path,
        mutant_fasta_path=mutant_fasta_path,
        dataset_type=dataset_type,
        window_size=window_size,
        strict_manifest=strict_manifest,
        progress=progress,
    )
