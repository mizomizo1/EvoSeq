from pathlib import Path
import tomllib

from .preprocess import preprocess_files, preprocess_folder
from .scoring import export_perbase_logprobs, score_pairs_file


def _none_if_blank(value):
    return None if value == "" else value


def load_config(path):
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def run_from_config(path):
    config = load_config(path)

    project = config.get("project", {})
    input_dir = project.get("input_dir", project.get("base_dir", "."))

    preprocess_config = config.get("preprocess", {})
    scoring_config = config.get("scoring", {})
    perbase_config = config.get("perbase", {})

    outputs = {}
    if preprocess_config.get("enabled", True):
        reference_fasta_path = _none_if_blank(preprocess_config.get("reference_fasta_path"))
        mutant_fasta_path = _none_if_blank(preprocess_config.get("mutant_fasta_path"))
        manifest_path = preprocess_config.get("manifest_path", "auto")

        if reference_fasta_path and mutant_fasta_path:
            evo_df, saved = preprocess_files(
                reference_fasta_path=reference_fasta_path,
                mutant_fasta_path=mutant_fasta_path,
                manifest_path=manifest_path,
                output_dir=_none_if_blank(
                    preprocess_config.get("output_dir", preprocess_config.get("out_dir"))
                ),
                strict_manifest=preprocess_config.get("strict_manifest", False),
                progress=preprocess_config.get("progress", True),
            )
        else:
            evo_df, saved = preprocess_folder(
                input_dir=input_dir,
                output_dir=_none_if_blank(
                    preprocess_config.get("output_dir", preprocess_config.get("out_dir"))
                ),
                manifest_path=manifest_path,
                reference_fasta_path=reference_fasta_path,
                mutant_fasta_path=mutant_fasta_path,
                dataset_type=preprocess_config.get("dataset_type", "auto"),
                window_size=preprocess_config.get("window_size"),
                strict_manifest=preprocess_config.get("strict_manifest", False),
                progress=preprocess_config.get("progress", True),
            )
        outputs["preprocess_df"] = evo_df
        outputs["preprocess_paths"] = saved

    if scoring_config.get("enabled", False):
        pairs_path = (
            _none_if_blank(scoring_config.get("pairs_path"))
            or outputs.get("preprocess_paths", {}).get("pairs")
        )
        if not pairs_path:
            raise ValueError(
                "scoring.enabled is true, but no pairs_path was provided and "
                "preprocessing did not produce a pair table."
            )

        result_df, paths = score_pairs_file(
            pairs_path=pairs_path,
            output_dir=_none_if_blank(
                scoring_config.get("output_dir", scoring_config.get("result_dir"))
            ),
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

    if perbase_config.get("enabled", False):
        output_path = export_perbase_logprobs(
            fasta_path=perbase_config["fasta_path"],
            output_path=_none_if_blank(perbase_config.get("output_path")),
            output_dir=_none_if_blank(perbase_config.get("output_dir")),
            model_name=perbase_config.get("model_name", "evo2_7b"),
            device=perbase_config.get("device", "cuda:0"),
            center=perbase_config.get("center", 4096),
            half_window=perbase_config.get("half_window", 320),
            local_path=_none_if_blank(perbase_config.get("local_path")),
            progress=perbase_config.get("progress", True),
        )
        outputs["perbase_path"] = output_path

    outputs["config_path"] = Path(path)
    return outputs
