import numpy as np
import pandas as pd


def read_manifest(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def summarize_manifest(df):
    return {
        "rows": len(df),
        "columns": df.shape[1],
        "unique_record_id": df["record_id"].nunique(),
        "unique_samples": df["sample"].nunique(),
        "unique_genes": df["gene"].nunique(),
        "unique_hgvs": df["hgvs"].nunique(),
        "unique_spdi": df["spdi"].nunique(),
    }


def aggregate_manifest(df):
    return (
        df.groupby("record_id")
        .agg(
            samples_joined=("sample", lambda x: ";".join(sorted(set(map(str, x))))),
            n_manifest_rows=("sample", "size"),
            gene_manifest=("gene", "first"),
            hgvs=("hgvs", "first"),
            annotation=("annotation", "first"),
            chrom=("chrom", "first"),
            pos1=("pos1", "first"),
            ref=("ref", "first"),
            alt=("alt", "first"),
            spdi=("spdi", "first"),
            note=(
                "note",
                lambda x: ";".join(sorted(set(map(str, x.dropna()))))
                if x.notna().any()
                else np.nan,
            ),
        )
        .reset_index()
    )