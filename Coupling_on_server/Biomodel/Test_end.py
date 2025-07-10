#!/usr/bin/env python
# coding: utf-8


import xarray as xr
import numpy as np
import os
import re
import sys
from glob import glob


diagfi_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
start_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
startfi_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"

output_filename = diagfi_path + "diagfi_global.nc"
if os.path.exists(output_filename):
    os.remove(output_filename)

# Concatenate diagfiNN.nc in diagfi_global.nc

files = glob(diagfi_path+"diagfi*.nc")
# files = [f for f in files if os.path.basename(f) != "diagfi.nc"]
# files = [f for f in files if os.path.basename(f) != "diagfi_global.nc"]

if len(files) > 0:
    files = sorted(files, key=lambda f: int(re.search(r'diagfi(\d+)\.nc', f).group(1)))

    datasets = []
    time_offset = 0

    all_vars = set()

    for file in files:
        try:
            ds = xr.open_dataset(file, decode_times=False)
            all_vars.update(ds.data_vars.keys())
            ds.close()
        except Exception as e:
            print(f"Erreur while reading {file}: {e}")

    for file in files:
        print(f"Reading {file}...")
        ds = xr.open_dataset(file, decode_times=False)

        if 'Time' not in ds:
            raise ValueError(f"'Time' variable is missing in {file}.")

        ds['Time'] = ds['Time'] + time_offset
        time_offset = float(ds['Time'].values[-1])

        # Ajouter les variables manquantes avec des NaN
        for var in all_vars:
            if var not in ds:
                shape = tuple(ds.dims[dim] for dim in ds[var].dims) if var in ds else (ds.dims['Time'],)
                ds[var] = (('Time',), np.full((ds.dims['Time'],), np.nan))

        datasets.append(ds)

    print("Files concatenation...")
    global_ds = xr.concat(datasets, dim='Time')

    # Sauvegarde
    output_filename = "diagfi_global.nc"
    global_ds.to_netcdf(diagfi_path+output_filename)
    print(f"Concatenated file saved as : {output_filename}")

print("Deleting diagfi*.nc...")

for f in files:
    if os.path.exists(f):
        os.remove(f)
        print(f"Deleted : {f}")

files2 = glob(start_path+"start*.nc")
files2 = [f for f in files2 if os.path.basename(f) != "start.nc"]
files3 = glob(startfi_path+"startfi*.nc")
files3 = [f for f in files3 if os.path.basename(f) != "startfi.nc"]

print("Deleting start*.nc...")

for f in files2:
    if os.path.exists(f):
        os.remove(f)
        print(f"Deleted : {f}")

print("Deleting startfi*.nc...")

for f in files3:
    if os.path.exists(f):
        os.remove(f)
        print(f"Deleted : {f}")

ds = xr.open_dataset(diagfi_path+"diagfi_global.nc", decode_times=False)

Tfus = 271.4

end_bool = False

latitudes = ds.latitude
longitudes = ds.longitude

ds_last_669 = ds.isel(Time=slice(-669, None)) # 1 martian year = 669 earth days

# Loop on physical points (lon,lat)
for i, lat in enumerate(latitudes):
    for j, lon in enumerate(longitudes):
        toce_value = ds_last_669['tslab1'].isel(latitude=i, longitude=j).mean(dim='Time').values

        if rnat < 0.5:
            if toce_value > Tfus: end_bool = True

ds.close()

if end_bool == True :
    print(" All the ocean is frozen --> end of the simulation ")
    sys.exit(1)

