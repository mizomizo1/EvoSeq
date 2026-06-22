from pathlib import Path
import tomllib

from .preprocess import preprocess_from_base_dir
from .scoring import score_evo2_pairs


def _none_if_blank(value):
    return None if value == "" else value


def load_config(path):
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def run_from_config(path):
    config = load_config(path)

    project = config.get("project", {})
    base_dir = project.get("base_dir", ".")

    preprocess_config = config.get("preprocess", {})
    scoring_config = config.get("scoring", {})

    outputs = {}
    if preprocess_config.get("enabled", True):
        evo_df, saved = preprocess_from_base_dir(
            base_dir=base_dir,
            out_dir=_none_if_blank(preprocess_config.get("out_dir")),
            manifest_path=preprocess_config.get("manifest_path", "auto"),
            reference_fasta_path=_none_if_blank(
                preprocess_config.get("reference_fasta_path")
            ),
            mutant_fasta_path=_none_if_blank(
                preprocess_config.get("mutant_fasta_path")
            ),
            dataset_type=preprocess_config.get("dataset_type", "auto"),
            window_size=preprocess_config.get("window_size"),
            strict_manifest=preprocess_config.get("strict_manifest", False),
            progress=preprocess_config.get("progress", True),
        )
        outputs["preprocess_df"] = evo_df
        outputs["preprocess_paths"] = saved

    if scoring_config.get("enabled", False):
        result_df, paths = score_evo2_pairs(
            base_dir=base_dir,
            pairs_path=_none_if_blank(scoring_config.get("pairs_path")),
            result_dir=_none_if_blank(scoring_config.get("result_dir")),
            manifest_path=scoring_config.get("manifest_path", "auto"),
            model_name=scoring_config.get("model_name", "evo2_7b"),
            device=scoring_config.get("device", "cuda:0"),
            local_path=_none_if_blank(scoring_config.get("local_path")),
            batch_size=scoring_config.get("batch_size", 8),
            force_reload=scoring_config.get("force_reload", False),
            require_recommended_gpu=scoring_config.get(
                "require_recommended_gpu",
                True,
            ),
            progress=scoring_config.get("progress", True),
        )
        outputs["scoring_df"] = result_df
        outputs["scoring_paths"] = paths

    outputs["config_path"] = Path(path)
    return outputs
