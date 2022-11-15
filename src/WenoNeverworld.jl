module WenoNeverworld

export neverworld_grid, weno_neverworld_simulation, standard_outputs!
export increase_simulation_Δt!, update_simulation_clock!

using CUDA
using Printf
using JLD2
using Oceananigans
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Units
using Oceananigans.Architectures: arch_array
using Oceananigans.Grids: peripheral_node, inactive_node, on_architecture
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, ImmersedBoundaryCondition, GridFittedBottom

include("weno_neverworld_utils.jl")
include("initial_conditions.jl")
include("weno_neverworld.jl")
include("weno_neverworld_outputs.jl")

end # module WenoNeverworld
