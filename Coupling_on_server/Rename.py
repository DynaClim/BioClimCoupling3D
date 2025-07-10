import sys
import os
import shutil

if len(sys.argv) != 2:
    print("Usage : python Rename.py <number of iterations>")
    sys.exit(1)

i = int(sys.argv[1]) + 0

# Renommage de diagfi.nc en diagfi{i}.nc
if os.path.exists("diagfi.nc"):
    shutil.move("diagfi.nc", f"diagfi{i}.nc")
    print(f"Renommé : diagfi.nc -> diagfi{i}.nc")
else:
    print("Fichier non trouvé : diagfi.nc")

# Renommage de start.nc en start{i}.nc
if os.path.exists("start.nc"):
    shutil.move("start.nc", f"start{i}.nc")
    print(f"Renommé : start.nc -> start{i}.nc")
else:
    print("Fichier non trouvé : start.nc")

if os.path.exists("startfi.nc"):
    shutil.move("startfi.nc", f"startfi{i}.nc")
    print(f"Renommé : startfi.nc -> startfi{i}.nc")
else:
    print("Fichier non trouvé : startfi.nc")

# Renommage de restart.nc en start.nc
if os.path.exists("restart.nc"):
    shutil.move("restart.nc", "start.nc")
    print("Renommé : restart.nc -> start.nc")
else:
    print("Fichier non trouvé : restart.nc")

if os.path.exists("restartfi.nc"):
    shutil.move("restartfi.nc", "startfi.nc")
    print("Renommé : restartfi.nc -> startfi.nc")
else:
    print("Fichier non trouvé : restartfi.nc")
