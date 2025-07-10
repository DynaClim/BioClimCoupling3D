import pickle

# Grid size (adjust to your case)
n_lat = 48
n_lon = 64
n_tasks = 50  # number of blocks to create

# All grid points as (i, j) pairs
all_points = [(i, j) for i in range(n_lat) for j in range(n_lon)]

# Split into blocks for parallel processing
blocks = [all_points[k::n_tasks] for k in range(n_tasks)]

# Save to file
with open("points_partition.pkl", "wb") as f:
    pickle.dump(blocks, f)

print(f"{n_tasks} blocks of points generated and saved.")
