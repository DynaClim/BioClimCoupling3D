# -*- coding: utf-8 -*-
"""
Creating corrk tables
Atmosphere of CO2, H2, CH4 and variable H2O
Warning: this script must be run as administrator (on Windows)
"""

import exo_k as xk
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import time
import sys
import os
import shutil

# Ensure MKS units are used
xk.Settings().set_mks(True)

# Paths

datapath = '/home/users/m/meyerfra/earlymars_test_couplage_global/datadir/corrk_data/'
datapath_LR = '/home/users/m/meyerfra/earlymars_test_couplage_global/LowRes/'
gases_path = '/home/users/m/meyerfra/earlymars_test_couplage_global/'

def parse_gases_def(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    n_gases = int(lines[0])
    gas_names = lines[1:1 + n_gases]
    ratios = [float(x) for x in lines[1 + n_gases:1 + 2 * n_gases]]

    return dict(zip(gas_names, ratios))

# Load gas ratios
gases = parse_gases_def(gases_path+"gases.def")

fG = gases.get('CH4', 0)
fH = gases.get('H2', 0)
fC = gases.get('CO2', 0)
fN = gases.get('N2', 0)

if fG < 0 or fH < 0 or fC < 0 or fN < 0 or fG+fC+fH+fN>1.0001 :
    print("Error : fG < 0 or fH < 0 or fC < 0 or fN < 0 or fG+fC+fH+fN>1.0001 ")
    sys.exit(1)

# Delete folder to avoid errors
shutil.move(datapath+"early_mars_h2_co2_ch4_n2_varh2o/Q.dat", gases_path+"Q.dat")
shutil.rmtree(datapath+"early_mars_h2_co2_ch4_n2_varh2o")



# Creating the LMDZ tables for fixed ratios of CH4, CO2, H2 and N2 and variable H2O

# IR range
molecs5 = ['CO2', 'CH4', 'H2', 'H2O']
suffix5 = 'IR_LR.ktable.exorem.SI.h5'
database5 = xk.Kdatabase(molecs5, suffix5, search_path=datapath_LR)
x_array = np.array([1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2])
vgas5 = {'H2O': 'background'}
bg_gas5 = {'CO2': fC, 'N2': 'background', 'CH4': fG, 'H2': fH}
mix_var_gas5_IR = database5.create_mix_ktable5d(vgas_comp=vgas5, bg_comp=bg_gas5, x_array=x_array)

# VI range
molecs6 = ['CO2', 'CH4', 'H2', 'H2O']
suffix6 = 'VI_LR.ktable.exorem.SI.h5'
database6 = xk.Kdatabase(molecs6, suffix6, search_path=datapath_LR)
vgas6 = {'H2O': 'background'}
bg_gas6 = {'CO2': fC, 'N2': 'background', 'CH4': fG, 'H2': fH}
mix_var_gas6_VI = database6.create_mix_ktable5d(vgas_comp=vgas6, bg_comp=bg_gas6, x_array=x_array)

# Save the tables
mix_var_gas5_IR.write_LMDZ(datapath + 'early_mars_h2_co2_ch4_n2_varh2o', band='IR')
mix_var_gas6_VI.write_LMDZ(datapath + 'early_mars_h2_co2_ch4_n2_varh2o', band='VI')

xk.finalize_LMDZ_dir(datapath + 'early_mars_h2_co2_ch4_n2_varh2o', mix_var_gas5_IR.Nw, mix_var_gas6_VI.Nw)

# Move Q.dat back
shutil.move(gases_path+"Q.dat",datapath+"early_mars_h2_co2_ch4_n2_varh2o/Q.dat")

shutil.rmtree(datapath+"early_mars_h2_co2_ch4_n2_varh2o/40x35")
os.makedirs(datapath+"early_mars_h2_co2_ch4_n2_varh2o/40x35", exist_ok=True)

shutil.copytree(datapath+"early_mars_h2_co2_ch4_n2_varh2o/VI35", datapath+"early_mars_h2_co2_ch4_n2_varh2o/40x35", dirs_exist_ok=True)
shutil.rmtree(datapath+"early_mars_h2_co2_ch4_n2_varh2o/VI35")

shutil.copytree(datapath+"early_mars_h2_co2_ch4_n2_varh2o/IR40", datapath+"early_mars_h2_co2_ch4_n2_varh2o/40x35", dirs_exist_ok=True)
shutil.rmtree(datapath+"early_mars_h2_co2_ch4_n2_varh2o/IR40")
