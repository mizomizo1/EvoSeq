from evoseq import run_from_config


outputs = run_from_config("evoseq.example.toml")
print(outputs.keys())
