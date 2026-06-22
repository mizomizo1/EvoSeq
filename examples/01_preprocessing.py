from evoseq.preprocess import preprocess_folder


evo_df, paths = preprocess_folder("test")

print(evo_df.shape)
print(paths)
