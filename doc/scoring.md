# Evo2 Scoring

`score_evo2_pairs` reads `evo2_input/evo2_pairs.tsv`, validates the sequences,
loads one Evo2 model, scores reference and mutant windows, and exports delta
scores.

```python
from evoseq.scoring import score_evo2_pairs

result_df, paths = score_evo2_pairs(
    base_dir="/content/drive/MyDrive/project/Model_7B_4096_POS",
    model_name="evo2_7b",
    batch_size=8,
)
```

## Colab GPU Checks

The scoring step prints the detected GPU name, CUDA version, package versions,
and GPU memory. By default, EvoSeq checks recommended memory for common Evo2
models:

- `evo2_7b`: at least 14 GB
- `evo2_20b`: at least 70 GB

Set `require_recommended_gpu=False` if you intentionally want to bypass this
guard.

## Model Cache

EvoSeq reuses the model already loaded in the Python process. If a different
model is requested later, it raises an error instead of silently loading another
large model into the same GPU runtime.

Restart the Colab runtime when switching model sizes.

## Outputs

- `evo2_variant_scores_unique.tsv`
- `evo2_variant_scores_manifest.tsv` when `data/manifest.tsv` exists
- `environment_info.tsv`
- `scoring_report.tsv`
