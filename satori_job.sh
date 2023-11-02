#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=1TB

# Upload modules
module purge all
module add spack
module add cuda/11.4
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

# MPI specific exports
export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

# Julia specific enviromental variables
export COMMON="/nobackup/users/lcbrock/"
export JULIA_DEPOT_PATH="${COMMON}/depot"
export JULIA_CUDA_MEMORY_POOL=none
export JULIA="${COMMON}/julia/julia"

# Profile specific variable
export JULIA_NVTX_CALLBACKS=gc

# Number of threads in SLURM mode
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

export RX=1
export RY=4

srun --mpi=pmi2 ./launch.sh $JULIA --check-bounds=no --project example/run_mpi.jl
