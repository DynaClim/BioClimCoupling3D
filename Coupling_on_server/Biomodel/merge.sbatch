#!/bin/bash
#SBATCH --job-name=making_NCmap
#SBATCH --output=logs/merge.out
#SBATCH --error=logs/merge.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Ajuste si tu fais du threading interne
#SBATCH --time=1:00:00
#SBATCH --partition=shared-cpu
#SBATCH --mem=100G

# Charger l'environnement Python si besoin (selon ton cluster)
# module load python/3.x
# source ~/myenv/bin/activate  # si environnement virtuel
cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

# Lancer la tâche Python pour l'ID correspondant
python merge_results.py
