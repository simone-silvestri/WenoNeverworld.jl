module WenoNeverworld

export NeverworldGrid, weno_neverworld_simulation, neverworld_simulation_seawater, standard_outputs!, checkpoint_outputs!
export increase_simulation_Δt!, update_simulation_clock!, run_simulation!

using CUDA
using KernelAbstractions: @kernel, @index
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

# Change Richardson number calculation for smoothness reasons
using Oceananigans.Operators
using Oceananigans.BuoyancyModels: ∂z_b
import Oceananigans.TurbulenceClosures: Riᶜᶜᶠ 
@inline Riᶜᶜᶠ(i, j, k, grid, velocities, tracers, buoyancy) = ℑxyᶜᶜᵃ(i, j, k, grid, ℑxyᶠᶠᵃ, local_Riᶜᶜᶠ, velocities, tracers, buoyancy)
@inline function local_Riᶜᶜᶠ(i, j, k, grid, velocities, tracers, buoyancy)
    ∂z_u² = ℑxᶜᵃᵃ(i, j, k, grid, ∂zᶠᶜᶠ, velocities.u)^2
    ∂z_v² = ℑyᵃᶜᵃ(i, j, k, grid, ∂zᶜᶠᶠ, velocities.v)^2
    N² = ∂z_b(i, j, k, grid, buoyancy, tracers)
    S² = ∂z_u² + ∂z_v²
    Ri = N² / S²

    # Clip N² and avoid NaN
    return ifelse(N² <= 0, zero(grid), Ri)
end

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

