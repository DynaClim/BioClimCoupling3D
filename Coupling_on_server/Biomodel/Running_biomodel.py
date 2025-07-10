#!/usr/bin/env python
# coding: utf-8

import os
import pickle
import numpy as np
import glob
import xarray as xr
from Physical_functions_GCM import *
from Constants_GCM import *
from Main_evolution_function_RK import *
from Bio_model import *
from Functions import *

# === PARAMETERS AND FUNCTIONS ===

def alphaC(T):
    """
    Computes the CO2 solubility depending on temperature

    Params:
     - T : float, temperature (K)

    Returns:
     - float, CO2 solubility (mol.L-1.bar-1)"""
    return np.exp(9345.17 / T - 167.8108 + 23.3585 * np.log(T) + (0.023517 - 2.3656e-4 * T + 4.7036e-7 * T**2) * 35.0)

def parse_gases_def(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    n_gases = int(lines[0])
    gas_names = lines[1:1 + n_gases]
    ratios = [float(x) for x in lines[1 + n_gases:1 + 2 * n_gases]]

    return dict(zip(gas_names, ratios))

# === TASK ID ===

task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))  # ou fallback à 0 si en local

# === PATHS ===
diagfi_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
gases_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
initialisation_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"
partition_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"

# === LOADIND POINTS BLOCK (allows parallelisation) ===
with open(partition_path+"points_partition.pkl", "rb") as f:
    blocks = pickle.load(f)

my_points = blocks[task_id]

# === LOADIND DATA ===
with xr.open_dataset(diagfi_path + "diagfi_global.nc", decode_times=False) as ds:
    ds = ds.load() 

times = ds["Time"]
duration = abs(times[-1] - times[-2])
n669 = int(np.ceil(669 / duration))
ds_last_669 = ds.isel(Time=slice(-n669, None)) # 1 martian year = 669 earth days

AREAS = ds["aire"].isel(Time=-1).values

with open(initialisation_path+"Initialisation_table.pkl", "rb") as f:
    Ini_data = pickle.load(f)
Ini_conc = Ini_data["initial_concentrations"]
tf = Ini_data["tf"]
natm_grid = Ini_data['natm_grid']


gases = parse_gases_def(gases_path+"gases.def")
fG = gases.get("CH4", 0)
fH = gases.get("H2", 0)
fC = gases.get("CO2", 0)
fN = gases.get("N2_", 0)

altitude = ds.altitude.values * 1e3
altitude = np.insert(altitude, 0, 0)
last_time_idx = -1

# === TRAITEMENT DES POINTS DU BLOC ===

results = {}

for i, j in my_points:

    if Ini_conc[i,j] == None:
        continue

    temp_profile = ds_last_669['temp'].isel(latitude=i, longitude=j).mean(dim='Time').values
    p_profile = ds_last_669['p'].isel(latitude=i, longitude=j).mean(dim='Time').values
    tsurf_value = ds_last_669['tsurf'].isel(latitude=i, longitude=j).mean(dim='Time').values
    toce_value = ds_last_669['tslab1'].isel(latitude=i, longitude=j).mean(dim='Time').values
    psurf_value = ds_last_669['ps'].isel(latitude=i, longitude=j).mean(dim='Time').values
    point_area = AREAS[i,j]

    temp_profile = np.insert(temp_profile, 0, tsurf_value)
    p_profile = np.insert(p_profile, 0, psurf_value)

    dico_ini = Ini_conc[i,j]
    NC = dico_ini['NC']
    cH = dico_ini['cH']
    cG = dico_ini['cG']
    cC = dico_ini['cC']
    cN = dico_ini['cN']

    Traits = ReturnTraits(toce_value, Par, 1)
    X0 = 2 * Traits[3]

    Bio, Medium,Atmo, t,flux,flux_times,Pression=system_evolution_RK(psurf_value,
                                                                     toce_value,
                                                                     550*point_area / fraction_ocean,
                                                                     1/11,
                                                                     tf, # simulation finish at tf
                                                                     point_area, 
                                                                     natm_grid[i,j],
                                                                     X0,
                                                                     CH4_frac_threshold=10, #10 > 1 => it is not the output condition
                                                                     focean=1, 
                                                                     Atm_compo=(fH, fC, fG,fN),
                                                                     concentration_ini=(cH,cC,cG,cN),
                                                                     NC_0=NC, 
                                                                     rtol=1e-7, 
                                                                     atol=1e-20,
                                                                     methode='LSODA',
                                                                     N=int(1e5),
                                                                     firststep=None,
                                                                     minstep=np.nan)
    results[(i, j)] = {
        "t": t,
        "Atmo": Atmo,
        "Medium": Medium,
        "Bio":Bio,
        "Pression":Pression,
        "timep":flux_times}

# === SAUVEGARDE DES RÉSULTATS ===

with open(f"results_{task_id}.pkl", "wb") as f:
    pickle.dump(results, f)

print(f"✅ Tâche {task_id} terminée. Résultats enregistrés.")





