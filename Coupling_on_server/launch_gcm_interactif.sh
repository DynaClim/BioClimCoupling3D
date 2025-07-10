#!/bin/bash

# Demande une session interactive avec 10 tâches sur le partition shared-cpu pendant 2h
salloc --ntasks=12 --partition=shared-cpu --time=08:00:00 <<'EOF'

module purge
module load GCC/11.2.0 OpenMPI/4.1.1 netCDF-Fortran/4.5.3 OpenBLAS/0.3.18

# Pour la compression bzip2
# export LD_LIBRARY_PATH=/opt/saltstack/salt/lib:$LD_LIBRARY_PATH

# Lancement du binaire
srun ./gcm_64x48x26_phystd_b40x35_para.e

EOF
