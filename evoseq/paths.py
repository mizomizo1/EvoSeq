from pathlib import Path


def common_parent(paths):
    resolved = [Path(path).expanduser().resolve().parent for path in paths if path]
    if not resolved:
        return Path.cwd()
    if len(resolved) == 1:
        return resolved[0]

    import os

    return Path(os.path.commonpath([str(path) for path in resolved]))


def default_output_dir(kind, *input_paths, base_dir=None):
    names = {
        "preprocess": "evoseq_preprocess_output",
        "scoring": "evoseq_scoring_output",
        "perbase": "evoseq_perbase_output",
    }
    dirname = names.get(kind, f"evoseq_{kind}_output")

    if base_dir:
        return Path(base_dir) / dirname

    return common_parent(input_paths) / dirname


def ensure_output_dir(path, fallback="/content/evoseq_output"):
    path = Path(path)
    try:
        path.mkdir(parents=True, exist_ok=True)
        test_path = path / ".write_test"
        test_path.write_text("ok")
        test_path.unlink(missing_ok=True)
        return path
    except OSError as exc:
        fallback_path = Path(fallback)
        print(f"Warning: cannot use output directory {path} ({exc}).")
        print(f"Using fallback output directory: {fallback_path}")
        fallback_path.mkdir(parents=True, exist_ok=True)
        return fallback_path
