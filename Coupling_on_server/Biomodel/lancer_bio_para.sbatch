#!/bin/bash
#SBATCH --job-name=bio_model
#SBATCH --output=logs/task_%A_%a.out
#SBATCH --error=logs/task_%A_%a.err
#SBATCH --array=0-49        # Lancer 50 tâches, ID de 0 à 49
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Ajuste si tu fais du threading interne
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --partition=shared-cpu

# Charger l'environnement Python si besoin (selon ton cluster)
# module load python/3.x
# source ~/myenv/bin/activate  # si environnement virtuel
cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel
# Lancer la tâche Python pour l'ID correspondant
python Running_biomodel.py
