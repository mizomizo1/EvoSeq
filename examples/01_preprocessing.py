from evoseq.preprocess import preprocess_from_base_dir


evo_df, paths = preprocess_from_base_dir(
    base_dir=".",
    dataset_type="auto",
    window_size=4096,
    manifest_path="auto",
)

print(evo_df.shape)
print(paths)
