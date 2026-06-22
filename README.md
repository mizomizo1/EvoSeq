# EvoSeq

EvoSeq is a small Colab-friendly toolkit for preparing paired reference/mutant
FASTA files and scoring variants with Evo2.

It is designed for the common research workflow where positive datasets have a
`manifest.tsv`, negative datasets may only have paired FASTA files, and the same
Evo2 model should stay loaded once per Colab runtime.

## Install

```bash
pip install -e .
```

For Evo2 scoring in Colab, install Evo2 and GPU dependencies in the runtime that
matches your model. The preprocessing step only needs the base dependencies.

## Quick Start: Preprocessing

If your project has `data/manifest.tsv` plus paired FASTA files under `data/`,
or only paired FASTA files for a negative set:

```python
from evoseq.preprocess import preprocess_from_base_dir

evo_df, paths = preprocess_from_base_dir(
    "/content/drive/MyDrive/my_project/Model_7B_4096_POS",
    dataset_type="auto",
    window_size=4096,
)
```

Outputs are written to `evo2_input/`:

- `evo2_pairs.tsv`: one row per variant with `ref_seq` and `mut_seq`
- `evo2_reference.fa`
- `evo2_mutant.fa`
- `evo2_all.fa`
- `preprocessing_report.tsv`

`manifest.tsv` is optional. When present, metadata are merged by `record_id`.
When absent, metadata are inferred from FASTA IDs when possible.

## Quick Start: Evo2 Scoring

```python
from evoseq.scoring import score_evo2_pairs

result_df, result_paths = score_evo2_pairs(
    base_dir="/content/drive/MyDrive/my_project/Model_7B_4096_POS",
    model_name="evo2_7b",
    batch_size=8,
)
```

Outputs are written to `evo2_results/`:

- `evo2_variant_scores_unique.tsv`
- `evo2_variant_scores_manifest.tsv` when a manifest is available
- `environment_info.tsv`
- `scoring_report.tsv`

Reference sequences are scored once per unique sequence and reused. This is
useful when many variants share the same reference window.

## Model Handling

EvoSeq caches the loaded Evo2 model inside the Python process:

```python
from evoseq.scoring import Evo2Scorer

scorer = Evo2Scorer(model_name="evo2_7b", device="cuda:0")
scores = scorer.score_sequences(["ACGTACGT"])
```

Calling another scoring function with the same model reuses it. Attempting to
load a different Evo2 model in the same runtime raises an explicit error by
default, because loading multiple large models often exhausts Colab GPU memory.
Restart the runtime when switching from 7B to 20B.

Common model names:

- `evo2_7b`
- `evo2_7b_base`
- `evo2_20b`

For local model weights:

```python
score_evo2_pairs(
    base_dir=".",
    model_name="evo2_20b",
    local_path="/content/drive/MyDrive/Models/evo2_20b.pt",
)
```

## TOML Config

Copy `evoseq.example.toml`, edit the paths/model, and run:

```python
from evoseq import run_from_config

outputs = run_from_config("evoseq.example.toml")
```

or:

```bash
evoseq-run evoseq.example.toml
```

## Reproducibility

EvoSeq writes small TSV reports for methods sections and reruns:

- input paths and output paths
- number of variants and unique reference sequences
- model name, batch size, device, elapsed time
- Python, PyTorch, CUDA, GPU, NumPy, pandas, Biopython, and Evo2 versions

These files are meant to be saved with each analysis directory.
