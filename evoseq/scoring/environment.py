from datetime import datetime
import platform
import sys


def collect_environment_info():
    info = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split(" ")[0],
        "platform": platform.platform(),
    }

    try:
        import torch

        info["torch"] = torch.__version__
        info["cuda_available"] = bool(torch.cuda.is_available())
        info["cuda_version"] = torch.version.cuda
        if torch.cuda.is_available():
            info["gpu_name"] = torch.cuda.get_device_name(0)
            info["gpu_total_memory_gb"] = round(
                torch.cuda.get_device_properties(0).total_memory / 1024**3,
                2,
            )
    except Exception as exc:
        info["torch_error"] = repr(exc)

    for package, module_name in [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("biopython", "Bio"),
        ("evo2", "evo2"),
    ]:
        try:
            module = __import__(module_name)
            info[package] = getattr(module, "__version__", "installed")
        except Exception:
            info[package] = "not installed"

    return info


def print_environment_info(info=None):
    info = info or collect_environment_info()
    print("=" * 60)
    print("Environment information")
    print("=" * 60)
    for key, value in info.items():
        print(f"{key:22}: {value}")
    print("=" * 60)


def write_environment_info(path, info=None):
    info = info or collect_environment_info()
    with open(path, "w") as fh:
        for key, value in info.items():
            fh.write(f"{key}\t{value}\n")
    return path
