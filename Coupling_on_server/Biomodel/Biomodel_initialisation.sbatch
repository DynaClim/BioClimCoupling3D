#!/bin/bash
#SBATCH --job-name=Bio_ini
#SBATCH --output=logs/Bio_ini.out
#SBATCH --error=logs/Bio_ini.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --partition=shared-cpu


cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

python Biomodel_initialisation.py --fGthresh $1
