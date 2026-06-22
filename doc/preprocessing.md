# Preprocessing

## Purpose

Create an Evo2-ready paired table and FASTA files from reference/mutant FASTA
files. A manifest is optional.

## Supported Inputs

- `data/manifest.tsv` with `record_id` for positive cohorts
- paired reference/mutant FASTA files
- FASTA IDs with either the project format
  `source|sample|gene|locus|variant` or simpler arbitrary IDs

## One Directory

```python
from evoseq.preprocess import preprocess_from_base_dir

evo_df, paths = preprocess_from_base_dir(
    "/content/drive/MyDrive/project/Model_7B_4096_POS",
    dataset_type="auto",
    window_size=4096,
)
```

This searches `data/` for manifest and FASTA files, then writes `evo2_input/`.

## Explicit Files

```python
from evoseq.preprocess import prepare_evo2_input

evo_df, paths = prepare_evo2_input(
    manifest_path=None,
    reference_fasta_path="data/evo2_neg_w16384_reference.fasta",
    mutant_fasta_path="data/evo2_neg_w16384_mutant.fasta",
    out_dir="evo2_input",
)
```

Set `manifest_path="data/manifest.tsv"` when a manifest is available.

## Outputs

- `evo2_pairs.tsv`
- `evo2_pair.tsv` for compatibility with older notebooks
- `evo2_reference.fa`
- `evo2_mutant.fa`
- `evo2_all.fa`
- `preprocessing_report.tsv`
