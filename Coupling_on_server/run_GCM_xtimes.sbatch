#!/bin/bash
#SBATCH --job-name=gcm_batch_loop
#SBATCH --output=output_loop_%j.log
#SBATCH --error=error_loop_%j.log
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --partition=shared-cpu

N=$1
if [ -z "$N" ]; then
  echo "Error: you must provide a number of iterations." >&2
  exit 1
fi

# Job ID of the previous step (initially empty)
prev_job=""

for ((i=0; i<N; i++)); do
  echo "------ Submitting iteration $i ------"

  # Submit the main executable job, with dependency if needed
  if [ -z "$prev_job" ]; then
    job1=$(sbatch --parsable run_executable.sbatch "$i")
  else
    job1=$(sbatch --parsable --dependency=afterok:$prev_job run_executable.sbatch "$i")
  fi

  if [ -z "$job1" ]; then
    echo "Error: failed to submit executable job at iteration $i" >&2
    exit 1
  fi

  # Submit the rename job, dependent on the executable job
  job2=$(sbatch --parsable --dependency=afterok:$job1 rename_xtimes.sbatch "$i")
  if [ -z "$job2" ]; then
    echo "Error: failed to submit rename job at iteration $i" >&2
    exit 1
  fi

  echo "Iteration $i: exec=$job1, rename=$job2"

  # The next iteration will depend on the renaming job of this one
  prev_job=$job2
done

echo "All jobs have been submitted with chained dependencies."
