#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH -t 12:00:00
#SBATCH --mem=90GB

module load cuda/11.4

JULIA="/home/sbishnu/Applications/julia/julia"

$JULIA --check-bounds=no --project simulations_quarter/simulation_seawater_quarter.jl
