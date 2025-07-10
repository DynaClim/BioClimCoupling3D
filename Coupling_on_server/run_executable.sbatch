#!/bin/bash
#SBATCH --job-name=2_gcm_run
#SBATCH --output=logs/gcm_%j.out
#SBATCH --error=logs/gcm_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=24
#SBATCH --partition=shared-cpu

module purge
module load GCC/11.2.0 OpenMPI/4.1.1 netCDF-Fortran/4.5.3 OpenBLAS/0.3.18

cd /home/users/m/meyerfra/earlymars_test_couplage_global

srun ./gcm_64x48x26_phystd_b40x35_para.e
