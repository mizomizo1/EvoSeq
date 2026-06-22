import numpy as np
import pandas as pd


def read_manifest(path):
    df = pd.read_csv(path, sep="\t")
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    if unnamed:
        df = df.drop(columns=unnamed)
    return df


def summarize_manifest(df):
    def nunique_or_zero(column):
        return df[column].nunique() if column in df.columns else 0

    return {
        "rows": len(df),
        "columns": df.shape[1],
        "unique_record_id": nunique_or_zero("record_id"),
        "unique_samples": nunique_or_zero("sample"),
        "unique_genes": nunique_or_zero("gene"),
        "unique_hgvs": nunique_or_zero("hgvs"),
        "unique_spdi": nunique_or_zero("spdi"),
    }


def aggregate_manifest(df, record_id_col="record_id"):
    if record_id_col not in df.columns:
        raise ValueError(f"Manifest is missing required column: {record_id_col}")

    grouped = df.groupby(record_id_col, dropna=False)
    out = grouped.size().rename("n_manifest_rows").reset_index()

    if record_id_col != "record_id":
        out = out.rename(columns={record_id_col: "record_id"})

    optional_first = [
        "gene",
        "hgvs",
        "annotation",
        "chrom",
        "pos1",
        "ref",
        "alt",
        "spdi",
    ]
    for column in optional_first:
        if column in df.columns:
            values = grouped[column].first().reset_index(drop=True)
            out[f"{column}_manifest" if column == "gene" else column] = values

    if "sample" in df.columns:
        out["samples_joined"] = grouped["sample"].agg(
            lambda x: ";".join(sorted(set(map(str, x.dropna()))))
        ).reset_index(drop=True)

    if "note" in df.columns:
        out["note"] = grouped["note"].agg(
            lambda x: ";".join(sorted(set(map(str, x.dropna()))))
            if x.notna().any()
            else np.nan
        ).reset_index(drop=True)

    return out
