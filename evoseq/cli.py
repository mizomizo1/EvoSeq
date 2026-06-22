import argparse

from .config import run_from_config


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run EvoSeq from a TOML config.")
    parser.add_argument("config", help="Path to an EvoSeq TOML config file.")
    args = parser.parse_args(argv)

    outputs = run_from_config(args.config)
    print("EvoSeq run completed.")
    for key, value in outputs.items():
        if key.endswith("_paths"):
            print(f"{key}:")
            for name, path in value.items():
                print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
