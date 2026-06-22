# Evo2 Scoring

`score_pairs_file` reads an EvoSeq pair table, validates the sequences, loads one
Evo2 model, scores reference and mutant windows, and exports delta scores.

```python
from evoseq.scoring import score_pairs_file

result_df, paths = score_pairs_file(
    pairs_path="test/evoseq_preprocess_output/evo2_pairs.tsv",
    model_name="evo2_7b",
    batch_size=8,
)
```

By default, results are written to an `evoseq_scoring_output/` directory next to
the pair table. Pass `output_dir` to choose a specific location.

## Required Runtime

Preprocessing does not require Evo2, PyTorch, or flash-attn. Real scoring does.
In Colab, install the runtime dependencies before scoring:

```bash
pip uninstall -y torchvision
pip install -q torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pip install -q flash-attn==2.8.0.post2 --no-build-isolation
pip install -q evo2
pip install -e .
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
