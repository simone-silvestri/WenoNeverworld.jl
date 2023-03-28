#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH -t 24:00:00
#SBATCH --mem=90GB

module load cuda/11.4

JULIA="/nobackup/users/hannahlu/julia/julia"

$JULIA --check-bounds=no --project run_high_resolution.jl 
