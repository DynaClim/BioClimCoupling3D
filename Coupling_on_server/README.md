# How to Launch the 3D Coupling on a Server?

This guide explains how to launch the 3D bio-climatic coupling using the GCM as the 3D climate model.

## 📑 Table of Contents

1. [Biomodel Directory](#1-biomodel-directory)  
2. [chemnetwork and datadir Directories](#2-chemnetwork-and-datadir-directories)  
3. [Required Files for the GCM](#3-required-files-for-the-gcm)  
4. [Files to Launch the GCM](#4-files-to-launch-the-gcm)  
5. [Creating Corrk Tables](#5-creating-corrk-tables)  
6. [Coupling](#6-coupling)  

## Initial Setup

First, install the GCM in your repository on the server.  
Refer to the GCM [documentation](https://lmdz-forge.lmd.jussieu.fr/mediawiki/Planets/index.php/Quick_Install_and_Run) for installation instructions.

Once installed, launch the following batch command in `LMDZ.COMMON` directory:

    ./makelmdz_fcm -arch BAO_YGG -parallel mpi  -p std -d 64x48x26 -b 32x36 -j 8 -full gcm

This will generate a file named: `gcm_64x48x26_phystd_b40x35_para.e`.

Then, copy this file into a new empty directory and make it executable:

    chmod +x gcm_64x48x26_phystd_b40x35_para.e

Next, copy all elements from this folder to the location of the executable file.  
At this point, the most difficult steps are completed. The only thing left is to update the paths in the necessary files.

> ## ⚠️ WARNING
>
> - The coupling has **never been fully run**, and the accuracy of the code is **not guaranteed**. However, one iteration of the coupling runs without issue.
> - The initial files representing the **prebiotic equilibrium state** are not yet generated. They must be created and added to the folder as:
>     - `diagfi_global.nc`
>     - `start.nc` (equivalent to `restart.nc` at the end of the GCM run)
>     - `startfi.nc` (equivalent to `restartfi.nc`)


## Diagram of the 3D coupling

![Server_3D_coupling.png](3D coupling diagram)

## 1. Biomodel Directory

This directory contains all the biological model codes that evolve the system.  
A `README.md` file is dedicated to it.

## 2. chemnetwork and datadir Directories

These directories **must not be modified**.  
They contain physical and chemical files essential for a successful GCM run.

⚠️ The `datadir` **must be used for coupling**, as it contains absorption tables tailored to the atmospheric composition (see `gases.def`).

## 3. Required Files for the GCM

### 3.1 Files not to be modified

- `Bands_64x48x26_24prc.dat`
- `num_run`
- `traceur.def`
- `z2sig.def`

### 3.2 `callphys.def`

This file allows modification of GCM physical parameters.  
Nothing to change. See [documentation](https://lmdz-forge.lmd.jussieu.fr/mediawiki/Planets/index.php/The_callphys.def_Input_File) for more details.  
Notable value: `Tfus = -2°C`, no photochemistry.

### 3.3 `run.def`

This file modifies the evolution parameters of the GCM.  
Nothing to change. See [documentation](https://lmdz-forge.lmd.jussieu.fr/mediawiki/Planets/index.php/The_run.def_Input_File).

### 3.4 `gases.def`

Defines the atmospheric composition.  
It evolves at each iteration due to biological changes.

Add a description

## 4. Files to Launch the GCM

### 4.1 `run_executable.sbatch`

This file launches the GCM in parallel. Update this line (path) :

    cd /home/users/m/meyerfra/earlymars_test_couplage_global

### 4.2 `Rename.py` and `rename_xtimes.sbatch`

These rename files after a GCM run to prepare the next one:

- `diagfi.nc` → `diagfi{iteration}.nc`
- `startfi.nc` → `startfi{iteration}.nc`
- `start.nc` → `start{iteration}.nc`
- `restart.nc` → `start.nc`
- `restartfi.nc` → `startfi.nc`

Run the batch with path update:

    cd /home/users/m/meyerfra/earlymars_test_couplage_global

### 4.3 `run_GCM_xtimes.sbatch`

This launches the GCM multiple times:

    sbatch run_GCM_xtimes.sbatch n

Where `n` is the number of iterations to run. 
It is not directly used in the coupling.
This is useful to generate the initial prebiotic state.

## 5. Creating Corrk Tables

### Creating_corrk_tables_gal.py

Creates corrk tables thanks to *Exok*.
Modify the following lines:

    datapath = '/home/users/m/meyerfra/earlymars_test_couplage_global/datadir/corrk_data/'
    datapath_LR = '/home/users/m/meyerfra/earlymars_test_couplage_global/LowRes/'
    gases_path = '/home/users/m/meyerfra/earlymars_test_couplage_global/'

### Corrktables.sbatch

Update these lines:

    cd /home/users/m/meyerfra/earlymars_test_couplage_global
    cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

## 6. Coupling

### 6.1 `list_fG.txt`

This list contains CH₄ thresholds to stop the biological model and re-run the GCM.  
Each value must be on a new line, no spaces, no comments, no blank lines.

A threshold corresponds to CH₄ threshol in one oceanic cell. To take into account atmospheric homogeneization, repeat each threshold 3 times. Example:

    1e-6
    1e-6
    1e-6
    3e-6
    3e-6
    3e-6
    ...

**⚠️ MUST BE EDITED BEFORE STARTING COUPLING!**

### 6.2 `Couplage_3D_global.sbatch`

Once parameters are fixed (e.g., `list_fG.txt`), run:

    sbatch Couplage_3D_global.sbatch

Update these variables in the script:

    LIST_FG_FILE="/home/users/m/meyerfra/earlymars_test_couplage_global/list_fG.txt"
    BASE_DIR="/home/users/m/meyerfra/earlymars_test_couplage_global"
    BIOMODEL_DIR="$BASE_DIR/Biomodel"
