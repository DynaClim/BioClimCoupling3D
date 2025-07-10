#!/bin/bash
#SBATCH --job-name=making_NCmap
#SBATCH --output=logs/merge.out
#SBATCH --error=logs/merge.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --time=1:00:00
#SBATCH --partition=shared-cpu
#SBATCH --mem=100G


cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

python merge_results.py
