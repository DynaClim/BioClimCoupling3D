#!/bin/bash
#SBATCH --job-name=CORRK
#SBATCH --output=logs/CORRK.out
#SBATCH --error=logs/CORRK.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --time=6:00:00
#SBATCH --mem=1G
#SBATCH --partition=shared-cpu

# Path
cd /home/users/m/meyerfra/earlymars_test_couplage_global

python Creating_corrk_tables_gal.py
ret=$?

# If the script fails (exit code ≠ 0), search for and cancel 'will never run' jobs
if [ $ret -ne 0 ]; then
    echo "Python script failed. Searching for and cancelling 'will never run' jobs..."

    cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel
    # Call the cleanup script
    ./cancel_unrunnable_jobs.sh
else
    echo "Python script completed successfully. No cleanup needed."
fi

