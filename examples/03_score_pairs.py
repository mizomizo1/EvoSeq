from evoseq.scoring import score_pairs_file


result_df, paths = score_pairs_file(
    pairs_path="test/evoseq_preprocess_output/evo2_pairs.tsv",
    model_name="evo2_7b",
    batch_size=8,
)

print(result_df[["id", "evo2_delta_score"]].head())
print(paths)
