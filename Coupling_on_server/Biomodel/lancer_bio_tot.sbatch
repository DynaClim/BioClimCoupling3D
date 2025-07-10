#!/bin/bash
#SBATCH --job-name=launch_bio_total
#SBATCH --output=logs/launch_total_bio.out
#SBATCH --error=logs/launch_total_bio.err
#SBATCH --time=12:00:00
#SBATCH --partition=shared-cpu

# ========== ARGUMENT PARSING ==========
# Use getopt to support --fG and --niter
PARSED_ARGS=$(getopt -o '' --long fG:,niter: -- "$@")
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to parse arguments."
    exit 1
fi
eval set -- "$PARSED_ARGS"

# Read options
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

# Check for required arguments
if [[ -z "$fG" || -z "$niter" ]]; then
    echo "Error: --fG and --niter arguments are required."
    exit 1
fi

# ========== PIPELINE LOGIC ==========
echo "===== Starting SLURM launch script ====="
echo "Arguments: fG=$fG | niter=$niter"

cd /home/users/m/meyerfra/earlymars_test_couplage_global/Biomodel

echo "Submitting generate_blocks.sbatch..."
job1=$(sbatch --parsable generate_blocks.sbatch)
echo "Job1 (generate_blocks) submitted with ID: $job1"

echo "Submitting Biomodel_initialisation.sbatch dependent on Job1..."
job2=$(sbatch --parsable --dependency=afterok:$job1 Biomodel_initialisation.sbatch $fG)
echo "Job2 (Biomodel_initialisation) submitted with ID: $job2"

echo "Submitting lancer_bio_para.sbatch dependent on Job2..."
job3=$(sbatch --parsable --dependency=afterok:$job2 lancer_bio_para.sbatch)
echo "Job3 (lancer_bio_para) submitted with ID: $job3"

echo "Submitting merge.sbatch dependent on Job3..."
job4=$(sbatch --parsable --dependency=afterok:$job3 merge.sbatch)
echo "Job4 (merge) submitted with ID: $job4"

echo "Submitting modify_files.sbatch dependent on Job4 with args fG=$fG and niter=$niter..."
job5=$(sbatch --parsable --dependency=afterok:$job4 modify_files.sbatch $niter)
echo "Job5 (modify_files) submitted with ID: $job5"

echo "===== End of SLURM launch script ====="
