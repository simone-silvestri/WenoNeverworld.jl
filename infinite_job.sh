#!/bin/bash

# Get the latest input file

while true; do
    # Check if the job is running
    if squeue -u $USER | grep -q "job_name"; then
        echo "Job is running"
    else
        # Restart the job with the latest input file
        latest_checkpoint=$(ls -t weno_thirtytwo* | head -n1)
        echo "rerunning simulation from output ${latest_checkpoint}"
        sed -i 's/init_file\s*=.*/init_file = "${latest_checkpoint}"/' run_mpi.jl
        sbatch -N1 satori_job.sh
    fi
    sleep 10
end
