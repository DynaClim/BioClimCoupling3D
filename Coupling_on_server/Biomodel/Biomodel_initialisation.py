#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np
import pickle
from IPython.utils import io
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
import sys
import argparse
import os

from Physical_functions_GCM import *
from Constants_GCM import *
from Bio_model import *
from Main_evolution_function_RK import *
from Functions import *

parser = argparse.ArgumentParser(description="f(CH4) threshold in the atmosphere")
parser.add_argument("--fGthresh", type=float, required=True, help="f(CH4) threshold in the atmosphere")
args = parser.parse_args()
fgt = args.fGthresh

print("fgt",fgt)

# Ocean fusion temperature

Tfus = 271.4    # K

# Defining paths

diagfi_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
gases_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/"
initialisation_path = "/home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel/"

# Extracting gases ratios in the atmosphere from gases.def

def parse_gases_def(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    n_gases = int(lines[0])
    gas_names = lines[1:1 + n_gases]
    ratios = [float(x) for x in lines[1 + n_gases:1 + 2 * n_gases]]

    return dict(zip(gas_names, ratios))

gases = parse_gases_def(gases_path+"gases.def")

if "CH4" in gases:
    fG0 = gases['CH4']
else: raise ValueError("CH4 not defined in gases.def")
if "H2" in gases:
    fH0 = gases['H2']
else: raise ValueError("H2 not defined in gases.def")
if "CO2" in gases:
    fC0 = gases['CO2']
else: raise ValueError("CO2 not defined in gases.def")
if "N2_" in gases:
    fN0 = gases['N2_']
else: raise ValueError("N2_ not defined in gases.def")

# ND interpolation for initial conditions

print("fH0,fC0,fG0,fN0",fH0,fC0,fG0,fN0)

df = pd.read_csv(initialisation_path+"Grid_values_ini.csv")

points = df[['P', 'T', 'fG']].values
NC_tab = df['NCeq'].values
cH_tab = df['cH'].values
cC_tab = df['cC'].values
cG_tab = df['cG'].values
cN_tab = df['cN'].values

interp_NCeq = LinearNDInterpolator(points, NC_tab)
interp_cH = LinearNDInterpolator(points, cH_tab)
interp_cC = LinearNDInterpolator(points, cC_tab)
interp_cG = LinearNDInterpolator(points, cG_tab)
interp_cN = LinearNDInterpolator(points, cN_tab)

def interp_conc(P, T, fG):
    pt = np.array([[P, T, fG]])
    return {
        'NC': interp_NCeq(pt)[0],
        'cH': interp_cH(pt)[0],
        'cC': interp_cC(pt)[0],
        'cG': interp_cG(pt)[0],
        'cN': interp_cN(pt)[0]
    }

# Loading diagi.nc
ds = xr.open_dataset(diagfi_path+"diagfi_global.nc", decode_times=False)

times = ds["Time"].values

last_time_idx = -1

latitudes = ds.latitude
longitudes = ds.longitude

altitude = ((ds.altitude).values) * 1e3
altitude = np.insert(altitude, 0, 0)

initial_concentrations = np.empty((len(latitudes), len(longitudes)), dtype=object)
natm_grid = np.zeros((len(latitudes), len(longitudes)))

duration = abs(times[-1] - times[-2])
n669 = int(np.ceil(669 / duration))
ds_last_669 = ds.isel(Time=slice(-n669, None)) # 1 martian year = 669 earth days

maxi = 0

rn = ds["rnat"].isel(Time=-1)
AREAS = ds["aire"].isel(Time=-1)

# Loop on physical points (lon,lat)
for i, lat in enumerate(latitudes):
    for j, lon in enumerate(longitudes):
        
        temp_profile = ds_last_669['temp'].isel(latitude=i, longitude=j).mean(dim='Time').values
        p_profile = ds_last_669['p'].isel(latitude=i, longitude=j).mean(dim='Time').values
        tsurf_value = ds_last_669['tsurf'].isel(latitude=i, longitude=j).mean(dim='Time').values
        toce_value = ds_last_669['tslab1'].isel(latitude=i, longitude=j).mean(dim='Time').values
        psurf_value = ds_last_669['ps'].isel(latitude=i, longitude=j).mean(dim='Time').values
        point_area = AREAS[i,j].values
        temp_profile = np.insert(temp_profile, 0, tsurf_value)
        p_profile = np.insert(p_profile, 0, psurf_value)
        rnat = rn[i,j]
        if rnat < 0.5:
            if toce_value <= Tfus: 
                initial_concentrations[i,j] = None
            else :
                initial_concentrations[i,j] = interp_conc(psurf_value, toce_value, fG0)
        else : initial_concentrations[i,j] = None

        npoint = atm_matter_quantity(p_profile, temp_profile,altitude, point_area, g=3.71)
        natm_grid[i,j] = npoint

        if point_area*interp_conc(psurf_value, toce_value, fG0)["NC"] > maxi:
            maxi = point_area*interp_conc(psurf_value, toce_value, fG0)["NC"]
            i_max,j_max = i,j

print("maxi=",maxi)
# Calculating tf for the physical point with the largest cells number
temp_profile_max = ds_last_669['temp'].isel(latitude=i_max, longitude=j_max).mean(dim='Time').values
p_profile_max = ds_last_669['p'].isel(latitude=i_max, longitude=j_max).mean(dim='Time').values
tsurf_value_max = ds_last_669['tsurf'].isel(latitude=i_max, longitude=j_max).mean(dim='Time').values
toce_value_max = ds_last_669['tslab1'].isel(latitude=i_max, longitude=j_max).mean(dim='Time').values
psurf_value_max = ds_last_669['ps'].isel(latitude=i_max, longitude=j_max).mean(dim='Time').values
point_area_max = AREAS[i_max,j_max].values

temp_profile_max = np.insert(temp_profile_max, 0, tsurf_value_max)
p_profile_max = np.insert(p_profile_max, 0, psurf_value_max)

print("toce,tsurf,psurf =",toce_value_max,tsurf_value_max,psurf_value_max)

dico_ini = initial_concentrations[i_max,j_max]
NC = dico_ini['NC']
cH = dico_ini['cH']
cG = dico_ini['cG']
cC = dico_ini['cC']
cN = dico_ini['cN']

print("NC,cH,cG,cN,cC",NC,cH,cG,cN,cC)

Traits = ReturnTraits(toce_value_max, Par, 1)
X0 = 2 * Traits[3]

Bio, Medium,Atmo, t,flux,flux_times,Pression=system_evolution_RK(psurf_value_max,
                                                                 toce_value_max,
                                                                 550*point_area_max / fraction_ocean,
                                                                 1/11,
                                                                 1e6,
                                                                 point_area_max, 
                                                                 natm_grid[i_max,j_max],
                                                                 X0,
                                                                 CH4_frac_threshold=fgt,
                                                                 focean=1, 
                                                                 Atm_compo=(fH0, fC0, fG0,fN0),
                                                                 concentration_ini=(cH,cC,cG,cN),
                                                                 NC_0=NC, 
                                                                 rtol=1e-7, 
                                                                 atol=1e-20,
                                                                 methode='LSODA',
                                                                 N=int(1e8),
                                                                 firststep=None,
                                                                 minstep=np.nan)

tf = t[-1]
print('tf=',tf)
print('t=',t)

if tf >= 9.999e5: raise ValueError("tf can not be defined --> Reduce CH4 threshold")
    
ds.close()

with open(initialisation_path+"Initialisation_table.pkl", "wb") as f:
    pickle.dump({
        'initial_concentrations': initial_concentrations,
        'natm_grid': natm_grid,
        'latitude': latitudes,
        'longitude': longitudes,
        'tf':tf
    }, f)

