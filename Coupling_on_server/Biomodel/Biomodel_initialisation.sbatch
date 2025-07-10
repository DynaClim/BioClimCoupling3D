#!/bin/bash
#SBATCH --job-name=Bio_ini
#SBATCH --output=logs/Bio_ini.out
#SBATCH --error=logs/bio_ini.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=50G
#SBATCH --partition=shared-cpu

# Charger l'environnement Python si besoin (selon ton cluster)
# module load python/3.x
# source ~/myenv/bin/activate  # si environnement virtuel
cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel
# Lancer la tâche Python pour l'ID correspondant
python Biomodel_initialisation.py --fGthresh $1
