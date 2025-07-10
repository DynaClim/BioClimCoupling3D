import pickle

# Taille de la grille (à adapter à ton cas)
n_lat = 48
n_lon = 64
n_tasks = 50  # nombre de blocs à créer

# Tous les points (i, j)
all_points = [(i, j) for i in range(n_lat) for j in range(n_lon)]

# Diviser en blocs
blocks = [all_points[k::n_tasks] for k in range(n_tasks)]

# Sauvegarder
with open("points_partition.pkl", "wb") as f:
    pickle.dump(blocks, f)

print(f"{n_tasks} blocs de points générés.")
