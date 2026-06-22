_MODEL_CACHE = {}
_ACTIVE_MODEL_KEY = None


class Evo2ModelAlreadyLoadedError(RuntimeError):
    pass


def _import_evo2_class():
    try:
        from evo2 import Evo2

        return Evo2
    except ImportError as first_exc:
        if getattr(first_exc, "name", None) not in {None, "evo2"}:
            raise RuntimeError(
                "Evo2 could not be imported because one of its runtime dependencies "
                f"is missing: {first_exc.name}. In Colab, install the Evo2 runtime "
                "dependencies before scoring, including torch, flash-attn, and evo2."
            ) from first_exc

    try:
        from evo2.models import Evo2

        return Evo2
    except ImportError as second_exc:
        raise RuntimeError(
            "Evo2 is not installed or is not importable. Preprocessing does not need "
            "Evo2, but scoring does. In Colab, install torch, flash-attn, and evo2 "
            "before calling score_pairs_file or export_perbase_logprobs."
        ) from second_exc


def ensure_cuda_device(device="cuda:0", require_gpu=True, min_memory_gb=None):
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError(
            "PyTorch is required for Evo2 scoring. Install the evoseq[evo2] extra "
            "or install torch in your Colab runtime."
        ) from exc

    if device.startswith("cuda") and not torch.cuda.is_available():
        message = "CUDA GPU is not available. In Colab, enable Runtime > Change runtime type > GPU."
        if require_gpu:
            raise RuntimeError(message)
        print(f"Warning: {message}")
        return device

    if device.startswith("cuda"):
        index = int(device.split(":", 1)[1]) if ":" in device else 0
        name = torch.cuda.get_device_name(index)
        props = torch.cuda.get_device_properties(index)
        total_gb = props.total_memory / 1024**3
        print(f"GPU detected: {name} ({total_gb:.1f} GB)")
        if min_memory_gb and total_gb < min_memory_gb:
            raise RuntimeError(
                f"GPU memory is {total_gb:.1f} GB, but {min_memory_gb:.1f} GB is requested."
            )

    return device


def load_evo2_model(
    model_name="evo2_7b",
    device="cuda:0",
    local_path=None,
    force_reload=False,
    prevent_cross_model_reload=True,
):
    global _ACTIVE_MODEL_KEY

    ensure_cuda_device(device=device, require_gpu=device.startswith("cuda"))

    key = (model_name, str(local_path) if local_path else None, device)
    if key in _MODEL_CACHE and not force_reload:
        print(f"Reusing loaded Evo2 model: {model_name}")
        return _MODEL_CACHE[key], device

    if _ACTIVE_MODEL_KEY is not None and _ACTIVE_MODEL_KEY != key and prevent_cross_model_reload:
        active_name, active_path, active_device = _ACTIVE_MODEL_KEY
        raise Evo2ModelAlreadyLoadedError(
            "An Evo2 model is already loaded in this Python process "
            f"({active_name}, local_path={active_path}, device={active_device}). "
            "Restart the runtime or pass force_reload=True after freeing GPU memory."
        )

    Evo2 = _import_evo2_class()

    print(f"Loading Evo2 model: {model_name}")
    if local_path:
        model = Evo2(model_name, local_path=str(local_path))
    else:
        model = Evo2(model_name)

    _MODEL_CACHE[key] = model
    _ACTIVE_MODEL_KEY = key
    return model, device


def tokenize_sequence(model, seq, device):
    import torch

    token_ids = model.tokenizer.tokenize(seq)
    return torch.tensor(token_ids, dtype=torch.int, device=device).unsqueeze(0)
