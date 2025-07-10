#!/bin/bash
#SBATCH --job-name=generate_blocks
#SBATCH --output=logs/generate_blocks.out
#SBATCH --error=logs/generate_blocks.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Ajuste si tu fais du threading interne
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --partition=shared-cpu

# Charger l'environnement Python si besoin (selon ton cluster)
# module load python/3.x
# source ~/myenv/bin/activate  # si environnement virtuel
cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel
# Lancer la tâche Python pour l'ID correspondant
python generate_points_blocks.py
