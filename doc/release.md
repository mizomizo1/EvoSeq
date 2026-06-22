# Release Guide

## Try Locally

From the repository root:

```bash
python -m pip install -e .
```

Then:

```python
from evoseq.preprocess import preprocess_folder

evo_df, paths = preprocess_folder("test")
print(paths)
```

For Evo2 scoring in Colab:

```bash
python -m pip install -e ".[evo2]"
```

## Build a Package

```bash
python -m pip install build
python -m build
```

This creates:

- `dist/evoseq-0.1.0.tar.gz`
- `dist/evoseq-0.1.0-py3-none-any.whl`

## GitHub Releases

GitHub Releases are a good first public distribution route.

1. Commit the release-ready code.
2. Create a tag, for example `v0.1.0`.
3. Build the package with `python -m build`.
4. Upload the wheel and source tarball from `dist/` to the GitHub Release.

Users can install from GitHub directly:

```bash
python -m pip install "git+https://github.com/hideakimizoue/EvoSeq.git@v0.1.0"
```

## PyPI

PyPI is better when you want users to install with `pip install evoseq`.

```bash
python -m pip install build twine
python -m build
python -m twine upload dist/*
```

Use TestPyPI first if you want a dry run.
