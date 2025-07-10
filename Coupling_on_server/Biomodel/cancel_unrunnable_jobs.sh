#!/bin/bash

# Optional: wait briefly to allow "will never run" jobs to appear
sleep 10

# Get all user jobs currently in the "DEPENDENCY" state
jobs_to_cancel=$(sacct --format=JobID,State --noheader | awk '$2=="DEPENDENCY" {print $1}')

if [ -z "$jobs_to_cancel" ]; then
    echo "No jobs to cancel."
else
    echo "Cancelling the following jobs:"
    echo "$jobs_to_cancel"
    scancel $jobs_to_cancel
fi
