module WenoNeverworld

export NeverworldGrid, weno_neverworld_simulation, standard_outputs!
export increase_simulation_Δt!, update_simulation_clock!, run_simulation!

using CUDA
using Printf
using JLD2
using Oceananigans
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Grids
using Oceananigans.Architectures: arch_array, architecture
using Oceananigans.Grids: on_architecture
using Oceananigans.ImmersedBoundaries

include("neverworld_bathymetry.jl")
include("neverworld_grid.jl")
include("neverworld_initial_and_boundary_conditions.jl")
include("weno_neverworld_utils.jl")
include("horizontal_visc.jl")
include("weno_neverworld.jl")
include("weno_neverworld_outputs.jl")

include("Diagnostics/Diagnostics.jl")

using .Diagnostics

end # module WenoNeverworld
