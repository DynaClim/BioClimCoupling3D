#!/bin/bash
#SBATCH --job-name=modify_files
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Ajuste si tu fais du threading interne
#SBATCH --time=1:00:00
#SBATCH --partition=shared-cpu
#SBATCH --mem=50G

cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel
exec > logs/modify_files_$1.out 2> logs/modify_files_$1.err

python Modifying_starting_files.py $1
