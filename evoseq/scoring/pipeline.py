from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .evo2_model import load_evo2_model
from .perbase import per_token_logprob, extract_center_window
from .export import read_fasta, write_perbase_logprobs_tsv
from .environment import (
    collect_environment_info,
    print_environment_info,
    write_environment_info,
)


RECOMMENDED_GPU_MEMORY_GB = {
    "evo2_7b": 14,
    "evo2_7b_base": 14,
    "evo2_20b": 70,
}


def _progress(iterable, enabled=True, **kwargs):
    if not enabled:
        return iterable
    try:
        from tqdm.auto import tqdm

        return tqdm(iterable, **kwargs)
    except Exception:
        return iterable


def _classify_variant(value):
    if isinstance(value, str) and ">" in value:
        ref, alt = value.split(">", 1)
        if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
            return "SNV"
    return "INDEL"


def validate_evo2_pairs(evo_df):
    required_columns = ["id", "ref_seq", "mut_seq"]
    missing = [column for column in required_columns if column not in evo_df.columns]
    if missing:
        raise ValueError(f"Missing required columns in Evo2 pair table: {missing}")

    evo_df = evo_df.copy()
    evo_df["ref_seq"] = evo_df["ref_seq"].astype(str).str.upper()
    evo_df["mut_seq"] = evo_df["mut_seq"].astype(str).str.upper()
    evo_df["ref_len"] = evo_df["ref_seq"].str.len()
    evo_df["mut_len"] = evo_df["mut_seq"].str.len()

    if "variant" not in evo_df.columns:
        evo_df["variant"] = evo_df["id"]
    if "gene" not in evo_df.columns:
        evo_df["gene"] = "NA"
    if "variant_type" not in evo_df.columns:
        evo_df["variant_type"] = evo_df["variant"].apply(_classify_variant)

    bad_ref = evo_df["ref_seq"].str.contains("[^ACGTN]", regex=True)
    bad_mut = evo_df["mut_seq"].str.contains("[^ACGTN]", regex=True)

    print("=" * 60)
    print("Evo2 input validation")
    print("=" * 60)
    print(f"Rows                    : {len(evo_df)}")
    print(f"Unique IDs              : {evo_df['id'].nunique()}")
    print(f"Duplicated IDs          : {evo_df['id'].duplicated().sum()}")
    print(f"Invalid ref sequences   : {bad_ref.sum()}")
    print(f"Invalid mutant sequences: {bad_mut.sum()}")
    print("Variant types:")
    print(evo_df["variant_type"].value_counts(dropna=False))
    print("Length pairs:")
    print(evo_df[["ref_len", "mut_len"]].value_counts().sort_index())

    if evo_df["id"].duplicated().any():
        raise ValueError("Duplicated id found. Evo2 scoring expects unique variants.")
    if bad_ref.any() or bad_mut.any():
        raise ValueError("Invalid sequence characters found.")

    return evo_df


def score_sequences_in_batches(model, sequences, batch_size=8, label="sequences", progress=True):
    scores = []
    n = len(sequences)
    ranges = range(0, n, batch_size)

    for start_idx in _progress(ranges, enabled=progress, desc=f"Scoring {label}"):
        end_idx = min(start_idx + batch_size, n)
        batch = sequences[start_idx:end_idx]
        print(f"Scoring {label}: {start_idx + 1}-{end_idx} / {n}")
        scores.extend(model.score_sequences(batch))

    return np.array(scores, dtype=float)


class Evo2Scorer:
    def __init__(
        self,
        model_name="evo2_7b",
        device="cuda:0",
        local_path=None,
        force_reload=False,
    ):
        self.model_name = model_name
        self.device = device
        self.local_path = local_path
        self.model, self.device = load_evo2_model(
            model_name=model_name,
            device=device,
            local_path=local_path,
            force_reload=force_reload,
        )
        self._executor = ThreadPoolExecutor(max_workers=1)

    def score_sequences(self, sequences, batch_size=8, label="sequences", progress=True):
        return score_sequences_in_batches(
            self.model,
            sequences,
            batch_size=batch_size,
            label=label,
            progress=progress,
        )

    def submit_score_sequences(self, sequences, batch_size=8, label="sequences"):
        return self._executor.submit(
            self.score_sequences,
            sequences,
            batch_size,
            label,
            False,
        )


def _resolve_pairs_path(base_dir=None, pairs_path=None):
    if pairs_path:
        return Path(pairs_path)
    if not base_dir:
        raise ValueError("Either base_dir or pairs_path is required.")
    return Path(base_dir) / "evo2_input" / "evo2_pairs.tsv"


def _resolve_result_dir(base_dir=None, result_dir=None):
    if result_dir:
        result_dir = Path(result_dir)
    elif base_dir:
        result_dir = Path(base_dir) / "evo2_results"
    else:
        result_dir = Path("evo2_results")
    result_dir.mkdir(parents=True, exist_ok=True)
    return result_dir


def _resolve_manifest_path(base_dir=None, manifest_path="auto"):
    if manifest_path in (None, False):
        return None
    if manifest_path != "auto":
        return Path(manifest_path)
    if not base_dir:
        return None
    candidate = Path(base_dir) / "data" / "manifest.tsv"
    return candidate if candidate.exists() else None


def score_evo2_pairs(
    base_dir=None,
    pairs_path=None,
    result_dir=None,
    manifest_path="auto",
    model_name="evo2_7b",
    device="cuda:0",
    local_path=None,
    batch_size=8,
    scorer=None,
    force_reload=False,
    require_recommended_gpu=True,
    progress=True,
):
    pairs_path = _resolve_pairs_path(base_dir=base_dir, pairs_path=pairs_path)
    result_dir = _resolve_result_dir(base_dir=base_dir, result_dir=result_dir)
    manifest_path = _resolve_manifest_path(base_dir=base_dir, manifest_path=manifest_path)

    min_memory_gb = (
        RECOMMENDED_GPU_MEMORY_GB.get(model_name)
        if require_recommended_gpu
        else None
    )
    if min_memory_gb:
        from .evo2_model import ensure_cuda_device

        ensure_cuda_device(device=device, min_memory_gb=min_memory_gb)

    env_info = collect_environment_info()
    print_environment_info(env_info)

    print(f"Evo2 pairs table: {pairs_path}")
    print(f"Result directory: {result_dir}")
    print(f"Model name      : {model_name}")
    print(f"Batch size      : {batch_size}")

    evo_df = pd.read_csv(pairs_path, sep="\t")
    evo_df = validate_evo2_pairs(evo_df)

    ref_unique = pd.Series(evo_df["ref_seq"].unique())
    ref_index = {seq: i for i, seq in enumerate(ref_unique)}
    ref_seq_indexes = evo_df["ref_seq"].map(ref_index).to_numpy()

    if scorer is None:
        scorer = Evo2Scorer(
            model_name=model_name,
            device=device,
            local_path=local_path,
            force_reload=force_reload,
        )

    start = time.time()
    print(f"Scoring {len(ref_unique)} unique reference sequences.")
    ref_scores = scorer.score_sequences(
        ref_unique.tolist(),
        batch_size=batch_size,
        label="reference",
        progress=progress,
    )
    ref_scores = np.asarray(ref_scores, dtype=float)

    print(f"Scoring {len(evo_df)} mutant sequences.")
    mut_scores = scorer.score_sequences(
        evo_df["mut_seq"].tolist(),
        batch_size=batch_size,
        label="mutant",
        progress=progress,
    )
    mut_scores = np.asarray(mut_scores, dtype=float)

    elapsed_sec = time.time() - start
    print(f"Scoring completed in {elapsed_sec:.1f} seconds.")

    result_unique = evo_df.copy()
    result_unique["evo2_ref_score"] = ref_scores[ref_seq_indexes]
    result_unique["evo2_mut_score"] = mut_scores
    result_unique["evo2_delta_score"] = (
        result_unique["evo2_mut_score"] - result_unique["evo2_ref_score"]
    )

    unique_result_path = result_dir / "evo2_variant_scores_unique.tsv"
    result_unique.to_csv(unique_result_path, sep="\t", index=False)
    print(f"Saved unique-variant scores: {unique_result_path}")

    manifest_result_path = None
    if manifest_path and manifest_path.exists():
        manifest_df = pd.read_csv(manifest_path, sep="\t")
        unnamed = [c for c in manifest_df.columns if str(c).startswith("Unnamed:")]
        if unnamed:
            manifest_df = manifest_df.drop(columns=unnamed)

        score_table = result_unique[
            [
                "id",
                "evo2_ref_score",
                "evo2_mut_score",
                "evo2_delta_score",
                "variant_type",
                "ref_len",
                "mut_len",
            ]
        ].rename(columns={"id": "record_id"})
        manifest_result = manifest_df.merge(score_table, on="record_id", how="left")
        manifest_result_path = result_dir / "evo2_variant_scores_manifest.tsv"
        manifest_result.to_csv(manifest_result_path, sep="\t", index=False)
        print(f"Saved manifest-level scores: {manifest_result_path}")
    else:
        print("Manifest was not found. Skipping manifest-level export.")

    env_path = result_dir / "environment_info.tsv"
    write_environment_info(env_path, env_info)

    report_path = result_dir / "scoring_report.tsv"
    report = {
        "pairs_path": str(pairs_path),
        "model_name": model_name,
        "local_path": str(local_path) if local_path else "",
        "device": device,
        "batch_size": batch_size,
        "n_variants": len(result_unique),
        "n_unique_reference_sequences": len(ref_unique),
        "elapsed_sec": round(elapsed_sec, 3),
        "mean_delta_score": result_unique["evo2_delta_score"].mean(),
        "median_delta_score": result_unique["evo2_delta_score"].median(),
        "std_delta_score": result_unique["evo2_delta_score"].std(),
        "unique_result_path": str(unique_result_path),
        "manifest_result_path": str(manifest_result_path or ""),
        "environment_info_path": str(env_path),
    }
    pd.DataFrame([report]).to_csv(report_path, sep="\t", index=False)
    print(f"Saved scoring report: {report_path}")

    return result_unique, {
        "unique_scores": unique_result_path,
        "manifest_scores": manifest_result_path,
        "environment_info": env_path,
        "scoring_report": report_path,
    }


def export_perbase_logprobs(
    fasta_path,
    output_path="perbase_logprobs.tsv",
    model_name="evo2_7b",
    device="cuda:0",
    center=4096,
    half_window=320,
    local_path=None,
    progress=True,
):
    model, device = load_evo2_model(
        model_name=model_name,
        device=device,
        local_path=local_path,
    )

    records = read_fasta(fasta_path)
    all_rows = []

    for tag, seq in _progress(records.items(), enabled=progress, desc="Per-base logprobs"):
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
