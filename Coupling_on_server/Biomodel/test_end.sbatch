#!/bin/bash
#SBATCH --job-name=TESTEND
#SBATCH --output=logs/TESTEND.out
#SBATCH --error=logs/TESTEND.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Ajuste si tu fais du threading interne
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --partition=shared-cpu

# Charger l'environnement Python si besoin (selon ton cluster)
# module load python/3.x
# source ~/myenv/bin/activate  # si environnement virtuel
cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

# Lancer la tâche Python pour l'ID correspondant
python Test_end.py
ret=$?

# Si le script échoue (code ≠ 0), alors chercher et annuler les jobs "will never run"
if [ $ret -ne 0 ]; then
    echo "Script Python a échoué. Recherche et suppression des jobs 'will never run'..."

    # Appel d’un script de nettoyage
    ./cancel_unrunnable_jobs.sh
else
    echo "Script Python a réussi. Aucun nettoyage nécessaire."
fi
