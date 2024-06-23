module WenoNeverworld

export NeverworldGrid, weno_neverworld_simulation, neverworld_simulation_seawater, standard_outputs!, checkpoint_outputs!
export WindStressBoundaryCondition, BuoyancyRelaxationBoundaryCondition
export increase_simulation_Δt!, update_simulation_clock!, run_simulation!
export years

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

const years = 365days

include("correct_oceananigans.jl")
include("Constants.jl")
include("Auxiliaries/Auxiliaries.jl")
include("NeverworldGrids/NeverworldGrids.jl")

include("NeverworldBoundaries/NeverworldBoundaries.jl")
include("Parameterizations/Parameterizations.jl")

include("weno_neverworld.jl")
include("weno_neverworld_outputs.jl")

include("Diagnostics/Diagnostics.jl")

using .Utils
using .NeverworldGrids
using .NeverworldBoundaries
using .Diagnostics

end # module WenoNeverworld

