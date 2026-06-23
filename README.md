# EvoSeq

<p align="center">
  <img src="doc/images/evoseq_banner.png" width="900">
</p>

EvoSeq is a small Colab-friendly toolkit for preparing paired reference/mutant
FASTA files and scoring variants with Evo2.

It is designed for workflows where positive datasets may include a `manifest.tsv`,
negative datasets may only have paired FASTA files, and the same Evo2 model should
stay loaded once per Colab runtime.

## Quick Start

### 1. Install

From PyPI:

```bash
pip install evoseq
```

For Evo2 scoring support:

```bash
pip install "evoseq[evo2]"
```

In Google Colab, Evo2 often needs runtime-specific GPU packages:

```bash
!pip uninstall -y torch torchvision torchaudio
!pip install -q torch==2.7.1 torchvision==0.22.1 torchaudio==2.7.1 --index-url https://download.pytorch.org/whl/cu128
!pip install -q flash-attn==2.8.0.post2 --no-build-isolation
!pip install -q evo2
!pip install -q evoseq
```

Preprocessing only needs the base `evoseq` install. Evo2 scoring requires a GPU
runtime with `torch`, `flash-attn`, and `evo2`.

### 2. Prepare paired FASTA files

Reference FASTA:

```fasta
>variant1
ACGTACGTACGT
```

Mutant FASTA:

```fasta
>variant1
ACGTTCGTACGT
```

The FASTA IDs must match between reference and mutant files.

### 3. Preprocess

```python
from evoseq.preprocess import preprocess_files

evo_df, paths = preprocess_files(
    reference_fasta_path="reference.fa",
    mutant_fasta_path="mutant.fa",
)
```

By default, outputs are written next to the input FASTA files in
`evoseq_preprocess_output/`.

Generated files:

* `evo2_pairs.tsv`: one row per variant with `ref_seq` and `mut_seq`
* `evo2_reference.fa`
* `evo2_mutant.fa`
* `evo2_all.fa`
* `preprocessing_report.tsv`

### 4. Score variants with Evo2

```python
from evoseq.scoring import score_pairs_file

result_df, result_paths = score_pairs_file(
    pairs_path=paths["pairs"],
    model_name="evo2_7b",
    batch_size=8,
)
```

By default, scoring outputs are written inside
`evoseq_preprocess_output/evoseq_scoring_output/`.

Generated files:

* `evo2_variant_scores_unique.tsv`
* `evo2_variant_scores_manifest.tsv` when a manifest is available
* `environment_info.tsv`
* `scoring_report.tsv`

Reference sequences are scored once per unique sequence and reused. This is useful
when many variants share the same reference window.

## Typical Workflow

```text
reference.fa + mutant.fa
            ↓
      preprocess
            ↓
      evo2_pairs.tsv
            ↓
        Evo2 scoring
            ↓
    variant score tables
```

Example:

```python
from evoseq.preprocess import preprocess_files
from evoseq.scoring import score_pairs_file

evo_df, paths = preprocess_files(
    reference_fasta_path="reference.fa",
    mutant_fasta_path="mutant.fa",
)

scores, outputs = score_pairs_file(
    pairs_path=paths["pairs"],
    model_name="evo2_7b",
)

print(scores.head())
```

Typical output structure:

```text
Typical output structure:

project/
├── reference.fa
├── mutant.fa
└── evoseq_preprocess_output/
    ├── evo2_pairs.tsv
    ├── evo2_reference.fa
    ├── evo2_mutant.fa
    ├── evo2_all.fa
    ├── preprocessing_report.tsv
    └── evoseq_scoring_output/
        ├── evo2_variant_scores_unique.tsv
        ├── evo2_variant_scores_manifest.tsv
        ├── scoring_report.tsv
        └── environment_info.tsv
```

The important files are:

* `evo2_pairs.tsv`: the main table passed to Evo2 scoring. It contains matched reference and mutant sequences.
* `evo2_variant_scores_unique.tsv`: the main Evo2 scoring result table.
* `evo2_variant_scores_manifest.tsv`: the score table merged with `manifest.tsv`, when a manifest is available.
* `preprocessing_report.tsv`: records what was generated during preprocessing.
* `scoring_report.tsv`: records model name, device, batch size, elapsed time, and output paths.
* `environment_info.tsv`: records package, CUDA, and GPU versions for reproducibility.


## Manifest Support

`manifest.tsv` is optional.

When present, metadata are merged by `record_id`:

```python
from evoseq.preprocess import preprocess_files

evo_df, paths = preprocess_files(
    reference_fasta_path="reference.fa",
    mutant_fasta_path="mutant.fa",
    manifest_path="manifest.tsv",
)
```

You can also let EvoSeq look for a manifest automatically:

```python
evo_df, paths = preprocess_files(
    reference_fasta_path="reference.fa",
    mutant_fasta_path="mutant.fa",
    manifest_path="auto",
)
```

When no manifest is provided, metadata are inferred from FASTA IDs when possible.


## Folder Discovery

If a folder contains paired FASTA files, EvoSeq can discover them:

```python
from evoseq.preprocess import preprocess_folder

evo_df, paths = preprocess_folder("test")
```

## Custom Output Directories

Preprocessing:

```python
evo_df, paths = preprocess_files(
    reference_fasta_path="reference.fa",
    mutant_fasta_path="mutant.fa",
    output_dir="outputs/preprocessing",
)
```

Scoring:

```python
from evoseq.scoring import score_pairs_file

result_df, result_paths = score_pairs_file(
    pairs_path="outputs/preprocessing/evo2_pairs.tsv",
    model_name="evo2_7b",
    output_dir="outputs/scoring",
)
```

## Per-Base Log-Probabilities

EvoSeq can export per-base Evo2 log-probabilities for aligned sequences:

```python
from evoseq.scoring import export_perbase_logprobs

path = export_perbase_logprobs(
    fasta_path="representative_perbase.fasta",
    model_name="evo2_7b",
    center=4096,
    half_window=320,
)
```

By default, this writes:

```text
evoseq_perbase_output/perbase_logprobs.tsv
```

The output can be used to visualize:

* reference vs mutant tracks
* Δ log-probability profiles
* local sequence effects
* long-range context effects

## Model Handling

EvoSeq caches the loaded Evo2 model inside the Python process:

```python
from evoseq.scoring import Evo2Scorer

scorer = Evo2Scorer(model_name="evo2_7b", device="cuda:0")
scores = scorer.score_sequences(["ACGTACGT"])
```

Calling another scoring function with the same model reuses it.

Attempting to load a different Evo2 model in the same runtime raises an explicit
error by default, because loading multiple large models often exhausts Colab GPU
memory. Restart the runtime when switching from 7B to 20B.

Common model names:

* `evo2_7b`
* `evo2_7b_base`
* `evo2_20b`

For local model weights:

```python
from evoseq.scoring import score_evo2_pairs

score_evo2_pairs(
    base_dir=".",
    model_name="evo2_20b",
    local_path="/content/drive/MyDrive/Models/evo2_20b.pt",
)
```

## TOML Config

Copy `evoseq.example.toml`, edit the input paths and model, then run:

```python
from evoseq import run_from_config

outputs = run_from_config("evoseq.example.toml")
```

or from the command line:

```bash
evoseq-run evoseq.example.toml
```

## Reproducibility

EvoSeq writes small TSV reports for methods sections and reruns.

Reports include:

* input paths and output paths
* number of variants and unique reference sequences
* model name, batch size, device, and elapsed time
* Python, PyTorch, CUDA, GPU, NumPy, pandas, Biopython, and Evo2 versions

Save these files alongside each analysis directory.

## Development

For local development from this repository:

```bash
git clone https://github.com/mizomizo1/EvoSeq.git
cd EvoSeq
pip install -e .
```

For local development with Evo2 extras:

```bash
pip install -e ".[evo2]"
```

Install a specific GitHub release directly:

```bash
pip install "git+https://github.com/mizomizo1/EvoSeq.git@v0.1.0"
```

Run tests without Evo2, torch, or flash-attn:

```bash
python -m unittest discover -s tests -v
```

These tests cover preprocessing, folder discovery, score-table export with a fake
scorer, and the missing Evo2 dependency message. Real Evo2 scoring requires a
Colab GPU runtime with `torch`, `flash-attn`, and `evo2` installed.
