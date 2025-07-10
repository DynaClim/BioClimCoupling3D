import sys
import os
import shutil

if len(sys.argv) != 2:
    print("Usage: python Rename.py <number of iterations>")
    sys.exit(1)

i = int(sys.argv[1]) + 0

# Rename diagfi.nc to diagfi{i}.nc
if os.path.exists("diagfi.nc"):
    shutil.move("diagfi.nc", f"diagfi{i}.nc")
    print(f"Renamed: diagfi.nc -> diagfi{i}.nc")
else:
    print("File not found: diagfi.nc")

# Rename start.nc to start{i}.nc
if os.path.exists("start.nc"):
    shutil.move("start.nc", f"start{i}.nc")
    print(f"Renamed: start.nc -> start{i}.nc")
else:
    print("File not found: start.nc")

# Rename startfi.nc to startfi{i}.nc
if os.path.exists("startfi.nc"):
    shutil.move("startfi.nc", f"startfi{i}.nc")
    print(f"Renamed: startfi.nc -> startfi{i}.nc")
else:
    print("File not found: startfi.nc")

# Rename restart.nc to start.nc
if os.path.exists("restart.nc"):
    shutil.move("restart.nc", "start.nc")
    print("Renamed: restart.nc -> start.nc")
else:
    print("File not found: restart.nc")

# Rename restartfi.nc to startfi.nc
if os.path.exists("restartfi.nc"):
    shutil.move("restartfi.nc", "startfi.nc")
    print("Renamed: restartfi.nc -> startfi.nc")
else:
    print("File not found: restartfi.nc")

