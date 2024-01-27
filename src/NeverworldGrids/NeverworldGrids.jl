module NeverworldGrids

using WenoNeverworld
using WenoNeverworld.Auxiliaries
using CUDA
using KernelAbstractions: @kernel, @index
using Printf
using JLD2
using Adapt
using Oceananigans
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.Architectures: arch_array, architecture
using Oceananigans.Grids: on_architecture
using Oceananigans.ImmersedBoundaries

export NeverworldGrid
export exponential_z_faces
export NeverWorldBathymetryParameters, neverworld_bathymetry

include("neverworld_bathymetry.jl")
include("neverworld_grid.jl")

end