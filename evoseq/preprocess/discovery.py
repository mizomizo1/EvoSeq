from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def _score_candidate(path, kind, dataset_type="auto", window_size=None):
    name = path.name.lower()
    score = 0

    if kind == "reference":
        score += 5 if "reference" in name else 0
        score += 3 if "_ref" in name or "ref_" in name else 0
    else:
        score += 5 if "mutant" in name else 0
        score += 3 if "_mut" in name or "mut_" in name else 0

    if dataset_type and dataset_type != "auto":
        aliases = {
            "positive": ["pos", "positive"],
            "negative": ["neg", "negative"],
        }.get(dataset_type, [dataset_type])
        score += 4 if any(alias in name for alias in aliases) else 0

    if window_size:
        score += 2 if str(window_size) in name else 0

    score -= 4 if "evo2_all" in name else 0
    score -= 3 if "output" in str(path).lower() else 0
    return score


def infer_dataset_type(base_dir):
    name = Path(base_dir).name.lower()
    if "neg" in name or "negative" in name:
        return "negative"
    if "pos" in name or "positive" in name:
        return "positive"
    return "auto"


def discover_manifest(base_dir, manifest_path="auto"):
    if manifest_path in (None, False):
        return None
    if manifest_path != "auto":
        return Path(manifest_path)

    base_dir = Path(base_dir)
    candidates = list((base_dir / "data").glob("manifest*.tsv"))
    candidates += list(base_dir.glob("manifest*.tsv"))
    return candidates[0] if candidates else None


def discover_fasta_pair(
    base_dir,
    reference_fasta_path=None,
    mutant_fasta_path=None,
    dataset_type="auto",
    window_size=None,
):
    if reference_fasta_path and mutant_fasta_path:
        return Path(reference_fasta_path), Path(mutant_fasta_path)

    base_dir = Path(base_dir)
    search_dirs = [base_dir / "data", base_dir]
    fasta_paths = []
    for search_dir in search_dirs:
        if search_dir.exists():
            fasta_paths.extend(
                p for p in search_dir.iterdir() if p.suffix.lower() in FASTA_SUFFIXES
            )

    if dataset_type == "auto":
        dataset_type = infer_dataset_type(base_dir)

    if not reference_fasta_path:
        refs = [
            p for p in fasta_paths
            if "reference" in p.name.lower() or "_ref" in p.name.lower()
        ]
        refs = sorted(
            refs,
            key=lambda p: _score_candidate(p, "reference", dataset_type, window_size),
            reverse=True,
        )
        reference_fasta_path = refs[0] if refs else None

    if not mutant_fasta_path:
        muts = [
            p for p in fasta_paths
            if "mutant" in p.name.lower() or "_mut" in p.name.lower()
        ]
        muts = sorted(
            muts,
            key=lambda p: _score_candidate(p, "mutant", dataset_type, window_size),
            reverse=True,
        )
        mutant_fasta_path = muts[0] if muts else None

    if not reference_fasta_path or not mutant_fasta_path:
        available = ", ".join(str(p) for p in sorted(fasta_paths)) or "none"
        raise FileNotFoundError(
            "Could not discover reference/mutant FASTA files. "
            "Pass reference_fasta_path and mutant_fasta_path explicitly. "
            f"Available FASTA files: {available}"
        )

    return Path(reference_fasta_path), Path(mutant_fasta_path)
