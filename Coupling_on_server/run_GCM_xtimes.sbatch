#!/bin/bash
#SBATCH --job-name=gcm_batch_loop
#SBATCH --output=output_loop_%j.log
#SBATCH --error=error_loop_%j.log
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --partition=shared-cpu

N=$1
if [ -z "$N" ]; then
  echo "Erreur : vous devez fournir un nombre d'itérations." >&2
  exit 1
fi

# Job ID de l'étape précédente (initialisé vide)
prev_job=""

for ((i=0; i<N; i++)); do
  echo "------ Soumission de l'itération $i ------"

  # Soumettre le job exécutable avec dépendance si nécessaire
  if [ -z "$prev_job" ]; then
    job1=$(sbatch --parsable run_executable.sbatch "$i")
  else
    job1=$(sbatch --parsable --dependency=afterok:$prev_job run_executable.sbatch "$i")
  fi

  if [ -z "$job1" ]; then
    echo "Erreur : soumission du job exécutable échouée à l'itération $i" >&2
    exit 1
  fi

  # Soumettre le job rename qui dépend du job exécutable
  job2=$(sbatch --parsable --dependency=afterok:$job1 rename_xtimes.sbatch "$i")
  if [ -z "$job2" ]; then
    echo "Erreur : soumission du job rename échouée à l'itération $i" >&2
    exit 1
  fi

  echo "Itération $i : exec=$job1, rename=$job2"

  # Pour la prochaine itération, on dépendra du renommage de cette itération
  prev_job=$job2
done

echo "Tous les jobs ont été soumis avec dépendances chaînées."
