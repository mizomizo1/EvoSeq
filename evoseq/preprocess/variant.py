import numpy as np
import pandas as pd


def compare_ref_mut(ref_records, mut_records, ids):
    rows = []

    for rid in sorted(ids):
        ref_seq = ref_records[rid]
        mut_seq = mut_records[rid]

        diffs = [
            i
            for i, (r, m) in enumerate(zip(ref_seq, mut_seq))
            if r != m
        ]

        rows.append(
            {
                "id": rid,
                "ref_len": len(ref_seq),
                "mut_len": len(mut_seq),
                "n_diffs": len(diffs),
                "first_diff_0based": diffs[0] if diffs else np.nan,
                "first_diff_1based": diffs[0] + 1 if diffs else np.nan,
                "ref_base": ref_seq[diffs[0]] if diffs else None,
                "mut_base": mut_seq[diffs[0]] if diffs else None,
            }
        )

    return pd.DataFrame(rows)


def parse_fasta_id(record_id):
    parts = record_id.split("|")

    if len(parts) < 5:
        return {
            "source": parts[0] if parts else "unknown",
            "sample_from_id": parts[1] if len(parts) > 1 else "NA",
            "gene": parts[2] if len(parts) > 2 else "NA",
            "locus": parts[3] if len(parts) > 3 else "NA",
            "variant": parts[4] if len(parts) > 4 else record_id,
        }

    return {
        "source": parts[0],
        "sample_from_id": parts[1],
        "gene": parts[2],
        "locus": parts[3],
        "variant": parts[4],
    }


def classify_variant(variant):
    if ">" not in variant:
        return "OTHER"

    ref, alt = variant.split(">", 1)

    if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
        return "SNV"

    return "INDEL"


def build_variant_table(pair_df):
    parsed = pair_df["id"].apply(parse_fasta_id).apply(pd.Series)
    variant_df = pd.concat([pair_df, parsed], axis=1)

    variant_df["variant_type"] = variant_df["variant"].apply(classify_variant)

    return variant_df
