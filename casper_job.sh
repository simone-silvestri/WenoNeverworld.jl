#!/bin/bash
#PBS -A UMIT0049
#PBS -N 4
#PBS -q casper@casper-pbs
#PBS -l walltime=01:00:00
#PBS -l select=2:ncpus=16:mpiprocs=16:ngpus=4:mem=300GB:mps=1
#PBS -l gpu_type=v100

# Use scratch for temporary files to avoid space limits in /tmp
export TMPDIR=/glade/scratch/ssilvest/temp
mkdir -p $TMPDIR

# Load modules to match compile-time environment
module purge
module load ncarenv nvhpc/22.5 cuda/11.4 openmpi/4.1.4

export CUDA_LAUNCH_BLOCKING=0 
export UCX_TLS=rc,sm,cuda_copy,cuda_ipc 
export OMPI_MCA_pml=ucx 
export OMPI_MCA_btl=self,vader,tcp,smcuda #openib 
export UCX_RNDV_SCHEME=get_zcopy 
export UCX_RNDV_THRESH=0 
export UCX_MAX_RNDV_RAILS=1 
export UCX_MEMTYPE_CACHE=n 

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

export JULIA="/glade/u/home/ssilvest/julia-1.9.3/bin/julia"

export RX=8
export RY=2

mpirun -n 16 ./launch.sh $JULIA --check-bounds=no --project scaling_experiments.jl 
