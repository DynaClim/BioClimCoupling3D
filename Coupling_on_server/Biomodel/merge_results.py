import pickle
import glob
import os

intermediate_dir = "merged_chunks"
os.makedirs(intermediate_dir, exist_ok=True)

files = sorted(glob.glob("results_*.pkl"))
chunk_size = 10
chunk_idx = 0

# First pass: process files in chunks to avoid memory overload
for i in range(0, len(files), chunk_size):
    chunk = files[i:i + chunk_size]
    combined = {}
    print(i, '/', len(files) - 1)
    for fname in chunk:
        try:
            with open(fname, "rb") as f:
                part = pickle.load(f)
                if isinstance(part, dict):
                    combined.update(part)
        except Exception as e:
            print(f"Error loading {fname}: {e}")

    with open(f"{intermediate_dir}/chunk_{chunk_idx}.pkl", "wb") as f_out:
        pickle.dump(combined, f_out)

    del combined
    chunk_idx += 1

# Second pass: merge all chunked files into one final result
final_result = {}
for chunk_file in sorted(glob.glob(f"{intermediate_dir}/chunk_*.pkl")):
    with open(chunk_file, "rb") as f:
        part = pickle.load(f)
        final_result.update(part)

with open("resultats_biomodel.pkl", "wb") as f:
    pickle.dump(final_result, f)

print("Merge completed successfully without memory overload.")
