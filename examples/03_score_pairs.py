from evoseq.scoring import score_evo2_pairs


result_df, paths = score_evo2_pairs(
    base_dir=".",
    model_name="evo2_7b",
    batch_size=8,
)

print(result_df[["id", "evo2_delta_score"]].head())
print(paths)
