#!/bin/bash
#SBATCH --job-name=rename_xtimes
#SBATCH --output=output_rename_sbatch_xtimes.log
#SBATCH --error=error_rename_sbatch_xtimes.log
#SBATCH --time=12:00:00
#SBATCH --partition=shared-cpu

i=$1

if [ -z "$i" ]; then
  echo "Error : no args . Utilisation : sbatch $0 <indice>"
  exit 1
fi

cd /home/users/m/meyerfra/earlymars_test_couplage_global

python Rename.py "$i"
