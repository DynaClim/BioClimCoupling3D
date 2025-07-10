#!/bin/bash
#SBATCH --job-name=bio_model
#SBATCH --output=logs/task_%A_%a.out
#SBATCH --error=logs/task_%A_%a.err
#SBATCH --array=0-49       
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --partition=shared-cpu

cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

python Running_biomodel.py
