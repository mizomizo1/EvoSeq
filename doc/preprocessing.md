# Preprocessing

## Purpose

Create an Evo2-ready paired table and FASTA files from reference/mutant FASTA
files. A manifest is optional.

## Supported Inputs

- `data/manifest.tsv` with `record_id` for positive cohorts
- paired reference/mutant FASTA files
- FASTA IDs with either the project format
  `source|sample|gene|locus|variant` or simpler arbitrary IDs

## Explicit Files

```python
from evoseq.preprocess import preprocess_files

evo_df, paths = preprocess_files(
    reference_fasta_path="test/evo2_reference.fasta",
    mutant_fasta_path="test/evo2_mutant.fasta",
    manifest_path="auto",
)
```

This writes `test/evoseq_preprocess_output/` unless `output_dir` is provided.

## Folder Discovery

```python
from evoseq.preprocess import preprocess_folder

evo_df, paths = preprocess_folder("test")
```

This searches the folder and its `data/` child for manifest and FASTA files.
Set `output_dir="outputs/preprocessing"` when you want a central output folder.

## Outputs

- `evo2_pairs.tsv`
- `evo2_pair.tsv` for compatibility with older notebooks
- `evo2_reference.fa`
- `evo2_mutant.fa`
- `evo2_all.fa`
- `preprocessing_report.tsv`
