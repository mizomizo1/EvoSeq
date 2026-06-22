# EvoSeq

EvoSeq is a small Colab-friendly toolkit for preparing paired reference/mutant
FASTA files and scoring variants with Evo2.

It is designed for the common research workflow where positive datasets have a
`manifest.tsv`, negative datasets may only have paired FASTA files, and the same
Evo2 model should stay loaded once per Colab runtime.

## Install

For local testing from this repository:

```bash
pip install -e .
```

For Evo2 scoring dependencies:

```bash
pip install -e ".[evo2]"
```

In Google Colab, Evo2 often needs a runtime-specific install. Use this before
scoring:

```bash
pip uninstall -y torchvision
pip install -q torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pip install -q flash-attn==2.8.0.post2 --no-build-isolation
pip install -q evo2
pip install -e .
```

After a GitHub Release is tagged, users can install a specific version directly:

```bash
pip install "git+https://github.com/mizomizo1/EvoSeq.git@v0.1.0"
```

For Evo2 scoring in Colab, install Evo2 and GPU dependencies in the runtime that
matches your model. The preprocessing step only needs the base dependencies.

## Debug / Test

Run the local workflow tests without Evo2, torch, or flash-attn:

```bash
python -m unittest discover -s tests -v
```

These tests cover preprocessing, folder discovery, score-table export with a
fake scorer, and the missing Evo2 dependency message. Real Evo2 scoring still
requires a Colab GPU runtime with `torch`, `flash-attn`, and `evo2` installed.

## Quick Start: Preprocessing Files

Put files anywhere, for example in `test/`, and pass the files directly:

```python
from evoseq.preprocess import preprocess_files

evo_df, paths = preprocess_files(
    reference_fasta_path="test/evo2_reference.fasta",
    mutant_fasta_path="test/evo2_mutant.fasta",
    manifest_path="auto",
)
```

By default, outputs are written next to the input files:
`test/evoseq_preprocess_output/`.

You can also be explicit:

```python
evo_df, paths = preprocess_files(
    reference_fasta_path="test/evo2_reference.fasta",
    mutant_fasta_path="test/evo2_mutant.fasta",
    output_dir="outputs/preprocessing",
)
```

Outputs include:

- `evo2_pairs.tsv`: one row per variant with `ref_seq` and `mut_seq`
- `evo2_reference.fa`
- `evo2_mutant.fa`
- `evo2_all.fa`
- `preprocessing_report.tsv`

`manifest.tsv` is optional. When present, metadata are merged by `record_id`.
When absent, metadata are inferred from FASTA IDs when possible.

## Quick Start: Preprocessing a Folder

If your folder contains paired FASTA files, EvoSeq can discover them:

```python
from evoseq.preprocess import preprocess_folder

evo_df, paths = preprocess_folder("test")
```

## Quick Start: Evo2 Scoring

```python
from evoseq.scoring import score_pairs_file

result_df, result_paths = score_pairs_file(
    pairs_path="test/evoseq_preprocess_output/evo2_pairs.tsv",
    model_name="evo2_7b",
    batch_size=8,
)
```

By default, outputs are written next to the pair table:
`test/evoseq_preprocess_output/evoseq_scoring_output/`.

Use `output_dir="outputs/scoring"` if you want a project-level result folder.

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

Copy `evoseq.example.toml`, edit the input paths/model, and run:

```python
from evoseq import run_from_config

outputs = run_from_config("evoseq.example.toml")
```

or:

```bash
evoseq-run evoseq.example.toml
```

## Per-Base Log-Probabilities

```python
from evoseq.scoring import export_perbase_logprobs

path = export_perbase_logprobs(
    fasta_path="test/representative_perbase.fasta",
    model_name="evo2_7b",
    center=4096,
    half_window=320,
)
```

By default, this writes `test/evoseq_perbase_output/perbase_logprobs.tsv`.

## Reproducibility

EvoSeq writes small TSV reports for methods sections and reruns:

- input paths and output paths
- number of variants and unique reference sequences
- model name, batch size, device, elapsed time
- Python, PyTorch, CUDA, GPU, NumPy, pandas, Biopython, and Evo2 versions

These files are meant to be saved with each analysis directory.
