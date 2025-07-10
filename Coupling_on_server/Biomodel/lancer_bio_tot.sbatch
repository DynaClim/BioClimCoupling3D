#!/bin/bash
#SBATCH --job-name=lancer_bio_tot
#SBATCH --output=logs/lancement_tot_bio.out
#SBATCH --error=logs/lancement_tot_bio.err
#SBATCH --time=12:00:00
#SBATCH --partition=shared-cpu

# ========== PARSING DES ARGUMENTS ==========
# Utilisation de getopt pour supporter --fG et --niter
PARSED_ARGS=$(getopt -o '' --long fG:,niter: -- "$@")
if [[ $? -ne 0 ]]; then
    echo "Erreur de parsing des arguments"
    exit 1
fi
eval set -- "$PARSED_ARGS"

# Lecture des options
while true; do
    case "$1" in
        --fG)
            fG="$2"
            shift 2
            ;;
        --niter)
            niter="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
    esac
done

# Vérification des arguments requis
if [[ -z "$fG" || -z "$niter" ]]; then
    echo "Erreur : les arguments --fG et --niter sont requis."
    exit 1
fi

# ========== LOGIQUE DU PIPELINE ==========
echo "===== Début du script de lancement SLURM ====="
echo "Arguments : fG=$fG | niter=$niter"

cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel


echo "Soumission de generate_blocks.sbatch..."
job1=$(sbatch --parsable generate_blocks.sbatch)
echo "Job1 (generate_blocks) soumis avec ID : $job1"

echo "Soumission de Biomodel_initialisation.sbatch dépendant de Job1..."
job2=$(sbatch --parsable --dependency=afterok:$job1 Biomodel_initialisation.sbatch $fG)
echo "Job2 (Biomodel_initialisation) soumis avec ID : $job2"

echo "Soumission de lancer_bio_para.sbatch dépendant de Job2..."
job3=$(sbatch --parsable --dependency=afterok:$job2 lancer_bio_para.sbatch)
echo "Job3 (lancer_bio_para) soumis avec ID : $job3"

echo "Soumission de merge.sbatch dépendant de Job3..."
job4=$(sbatch --parsable --dependency=afterok:$job3 merge.sbatch)
echo "Job4 (merge) soumis avec ID : $job4"

echo "Soumission de modify_files.sbatch dépendant de Job4 avec args fG=$fG niter=$niter..."
job5=$(sbatch --parsable --dependency=afterok:$job4 modify_files.sbatch $niter)
echo "Job5 (modify_files) soumis avec ID : $job5"

echo "===== Fin du script de lancement SLURM ====="
