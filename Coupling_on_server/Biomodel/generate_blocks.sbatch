#!/bin/bash
#SBATCH --job-name=generate_blocks
#SBATCH --output=logs/generate_blocks.out
#SBATCH --error=logs/generate_blocks.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --partition=shared-cpu


cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

python generate_points_blocks.py
