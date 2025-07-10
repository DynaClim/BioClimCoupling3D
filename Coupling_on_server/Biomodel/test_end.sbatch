#!/bin/bash
#SBATCH --job-name=TESTEND
#SBATCH --output=logs/TESTEND.out
#SBATCH --error=logs/TESTEND.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --partition=shared-cpu


cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

python Test_end.py
ret=$?


if [ $ret -ne 0 ]; then
    echo "Script Python failed. Searching and deleting jobs 'will never run'..."

    ./cancel_unrunnable_jobs.sh
else
    echo "Script Python successed. "
fi
