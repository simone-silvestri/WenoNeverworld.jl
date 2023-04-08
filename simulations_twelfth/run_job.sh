#!/bin/bash 
#SBATCH -N 1
#SBATCH -p sched_mit_raffaele_gpu 
#SBATCH --gres=gpu:1 
#SBATCH --ntasks-per-node=1
#SBATCH -t 24:00:00
##SBATCH --mem=80G
/home/sandre/Julia/julia-1.8.5/bin/julia --project --check-bounds=no twelfth.jl 