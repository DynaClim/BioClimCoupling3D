#!/bin/bash
#SBATCH --job-name=3Dcoupling
#SBATCH --output=3Dcoupling.out
#SBATCH --error=3Dcoupling.err
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-cpu

# ========== Paths ==========
LIST_FG_FILE="/home/users/m/meyerfra/earlymars_test_couplage_global/list_fG.txt"
BASE_DIR="/home/users/m/meyerfra/earlymars_test_couplage_global"
BIOMODEL_DIR="$BASE_DIR/Biomodel"

# ========== Read and display fG values ==========
# Read all non-empty lines into an array
mapfile -t fG_list < <(grep -v '^\s*$' "$LIST_FG_FILE")

# Display the list of fG values
echo "List of fG values:"
for fg in "${fG_list[@]}"; do
    echo " - $fg"
done
echo "Total number of fG values: ${#fG_list[@]}"

# ========== Initialize ==========
iteration_index=0
prev_chain_last_job=""

# ========== Loop over fG values ==========
for fG in "${fG_list[@]}"; do
    echo
    echo "=== Submitting job chain for fG=$fG (iteration $iteration_index) ==="

    # Optional dependency on the previous chain
    dep_flag=""
    if [ -n "$prev_chain_last_job" ]; then
        dep_flag="--dependency=afterok:$prev_chain_last_job"
        echo "Chaining after previous job ID: $prev_chain_last_job"
    fi

    # --------- BIOMODEL CHAIN ---------
    cd "$BIOMODEL_DIR" || exit 1

    echo "Submitting generate_blocks.sbatch..."
    job1=$(sbatch --parsable $dep_flag generate_blocks.sbatch)
    echo "  → Job ID: $job1"

    echo "Submitting Biomodel_initialisation.sbatch with fG=$fG..."
    job2=$(sbatch --parsable --dependency=afterok:$job1 Biomodel_initialisation.sbatch "$fG")
    echo "  → Job ID: $job2"

    echo "Submitting lancer_bio_para.sbatch..."
    job3=$(sbatch --parsable --dependency=afterok:$job2 lancer_bio_para.sbatch)
    echo "  → Job ID: $job3"

    echo "Submitting merge.sbatch..."
    job4=$(sbatch --parsable --dependency=afterok:$job3 merge.sbatch)
    echo "  → Job ID: $job4"

    echo "Submitting modify_files.sbatch with niter=$iteration_index..."
    job5=$(sbatch --parsable --dependency=afterok:$job4 modify_files.sbatch "$iteration_index")
    echo "  → Job ID: $job5"

    # --------- CORRKTABLES ---------
    echo "Submitting Corrktables.sbatch..."
    job6=$(sbatch --parsable --dependency=afterok:$job5 "$BASE_DIR/Corrktables.sbatch")
    echo "  → Job ID: $job6"

    # --------- GCM CHAIN (REPLACING run_GCM_xtimes) ---------
    num_gcm_iters=20
    prev_gcm_job=""

    cd "$BASE_DIR" || exit 1

    for ((gcm_i=0; gcm_i<num_gcm_iters; gcm_i++)); do
        echo "Submitting GCM iteration $gcm_i..."

        # Submit run_executable with dependency on previous GCM iteration or job6
        if [ -z "$prev_gcm_job" ]; then
            exec_job=$(sbatch --parsable --dependency=afterok:$job6 "$BASE_DIR/run_executable.sbatch" "$gcm_i")
        else
            exec_job=$(sbatch --parsable --dependency=afterok:$prev_gcm_job "$BASE_DIR/run_executable.sbatch" "$gcm_i")
        fi
        echo "  → run_executable.sbatch submitted: $exec_job"

        # Submit rename_xtimes depending on the executable job
        rename_job=$(sbatch --parsable --dependency=afterok:$exec_job "$BASE_DIR/rename_xtimes.sbatch" "$gcm_i")
        echo "  → rename_xtimes.sbatch submitted: $rename_job"

        # Update dependency for the next iteration
        prev_gcm_job="$rename_job"
    done

    # --------- FINAL TEST ---------
    echo "Submitting test_end.sbatch after all GCM iterations..."
    job8=$(sbatch --parsable --dependency=afterok:$prev_gcm_job "$BIOMODEL_DIR/test_end.sbatch")
    echo "  → Job ID: $job8"

    # Save last job of this chain to depend on in next iteration
    prev_chain_last_job="$job8"

    echo "Job chain for fG=$fG completed: $job1 → $job2 → $job3 → $job4 → $job5 → $job6 → ... → $rename_job → $job8"
    echo "==========================================================================="

    ((iteration_index++))
done
