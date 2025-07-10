#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import pickle
import numpy as np
import os
import shutil
import glob
import sys

# Paths
diagfi_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
initialisation_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"
biomodel_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"
gases_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
start_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
new_folder_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
results_xx_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"

n_iteration = int(sys.argv[1])
new_folder_name = os.path.join(new_folder_path, f"RESULTS_iteration{n_iteration}")
os.makedirs(new_folder_name, exist_ok=True)

# Utilitaires
def parse_gases_def(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    n_gases = int(lines[0])
    gas_names = lines[1:1 + n_gases]
    ratios = [float(x) for x in lines[1 + n_gases:1 + 2 * n_gases]]
    return dict(zip(gas_names, ratios))

gases = parse_gases_def(os.path.join(gases_path, "gases.def"))

previous_fG = gases.get("CH4", 0)
previous_fH = gases.get("H2", 0)
previous_fC = gases.get("CO2", 0)
previous_fN = gases.get("N2", 0)

# Coordonnees depuis diagfi.nc
ds = xr.open_dataset(os.path.join(diagfi_path, "diagfi_global.nc"), decode_times=False)
latitudes = ds.latitude.values
longitudes = ds.longitude.values
ds.load()
ds.close()

# Chargement des fichiers pickle
with open(os.path.join(biomodel_path, "resultats_biomodel.pkl"), "rb") as f:
    results = pickle.load(f)

with open(os.path.join(initialisation_path, "Initialisation_table.pkl"), "rb") as f:
    Ini_file = pickle.load(f)
    natm_map = Ini_file['natm_grid']

# Calcul des nouveaux facteurs
previous_ntotal = ntotal = ntotH = ntotC = ntotG = ntotN = 0

for i in range(len(latitudes)):
    for j in range(len(longitudes)):
        previous_ntotal += natm_map[i, j]
        if (i, j) in results:
            atmo = results[(i, j)]["Atmo"]
            last = atmo[:, -1]
            total = np.sum(last)
            ntotal += total
            ntotH += atmo[0, -1]
            ntotC += atmo[1, -1]
            ntotG += atmo[2, -1]
            ntotN += atmo[3, -1]
        else:
            ntotal += natm_map[i, j]
            ntotH += previous_fH * natm_map[i, j]
            ntotC += previous_fC * natm_map[i, j]
            ntotG += previous_fG * natm_map[i, j]
            ntotN += previous_fN * natm_map[i, j]

fH = ntotH / ntotal
fC = ntotC / ntotal
fG = ntotG / ntotal
fN = ntotN / ntotal

pressure_factor = ntotal / previous_ntotal

# Mise à jour de gases.def
def update_gas_concentrations(file_path, new_concentrations):
    with open(file_path, "r") as f:
        lines = f.readlines()
    gas_count = int(lines[1].strip())
    gas_names = [lines[i].strip() for i in range(2, 2 + gas_count)]
    ratio_start = 2 + gas_count
    updated_ratios = [f"{new_concentrations.get(gas, lines[i].strip())}\n" for i, gas in enumerate(gas_names, start=ratio_start)]
    updated_lines = lines[:ratio_start] + updated_ratios
    with open(file_path, "w") as f:
        f.writelines(updated_lines)

new_ratios = {"CO2": fC, "H2": fH, "CH4": fG, "N2_": fN}
shutil.copy2(os.path.join(gases_path, "gases.def"), os.path.join(new_folder_name, f'gases{n_iteration}.def'))
update_gas_concentrations(os.path.join(gases_path, "gases.def"), new_ratios)

# Mise à jour de start.nc
input_file = os.path.join(start_path, "start.nc")
ds2 = xr.open_dataset(input_file)
ds2.load()
ds2.close()

if "ps" in ds2 and ds2["ps"].size > 0:
    ds2["ps"] *= pressure_factor
    print("Variable 'ps' modifiée.")
else:
    print("Variable 'ps' absente ou vide.")

if "presnivs" in ds2 and ds2["presnivs"].size > 0:
    ds2["presnivs"] *= pressure_factor
    print("Variable 'presnivs' modifiée.")
else:
    print("Variable 'presnivs' absente ou vide.")

ds2.to_netcdf(input_file,mode='w')
print(f"Fichier modifié sauvegardé sous '{input_file}'.")
ds2.close()

# changing the names
shutil.copy(start_path + 'start.nc', new_folder_name + f'/start{n_iteration}.nc')
shutil.copy(start_path + 'startfi.nc', new_folder_name + f'/startfi{n_iteration}.nc')
os.rename(diagfi_path+'diagfi_global.nc', new_folder_name+f'/diagfi_global{n_iteration}.nc')
os.rename(biomodel_path+'resultats_biomodel.pkl', new_folder_name+f'/resultats_biomodel{n_iteration}.pkl')
os.rename(initialisation_path+'Initialisation_table.pkl', new_folder_name+f'/Initialisation_table{n_iteration}.pkl')

# Suppression sécurisée des fichiers results_*.pkl
fichiers = glob.glob(os.path.join(results_xx_path, "results_*.pkl"))
for fichier in fichiers:
    try:
        os.remove(fichier)
        print(f"Supprimé : {fichier}")
    except PermissionError:
        print(f"Permission refusée pour : {fichier}")
    except Exception as e:
        print(f"Erreur lors de la suppression de {fichier} : {e}")


