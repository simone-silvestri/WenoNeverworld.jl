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
using Oceananigans.Architectures: on_architecture, architecture
using Oceananigans.Grids: on_architecture
using Oceananigans.ImmersedBoundaries

export NeverworldGrid
export exponential_z_faces
export NeverworldBathymetry

include("neverworld_bathymetry.jl")
include("dino_grid_parameters.jl")
include("neverworld_grid.jl")

end