using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb
using Oceananigans.BoundaryConditions: ContinuousBoundaryFunction

import Oceananigans.Utils: getregion, _getregion

@inline getregion(cf::ContinuousBoundaryFunction{X, Y, Z, I}, i) where {X, Y, Z, I} =
    ContinuousBoundaryFunction{X, Y, Z, I}(cf.func,
                                           _getregion(cf.parameters, i),
                                           cf.field_dependencies,
                                           cf.field_dependencies_indices,
                                           cf.field_dependencies_interp)

@inline _getregion(cf::ContinuousBoundaryFunction{X, Y, Z, I}, i) where {X, Y, Z, I} =
    ContinuousBoundaryFunction{X, Y, Z, I}(cf.func,
                                           getregion(cf.parameters, i),
                                           cf.field_dependencies,
                                           cf.field_dependencies_indices,
                                           cf.field_dependencies_interp)

output_dir    = joinpath(@__DIR__, "../files_sixteen_new_bathy")
@show output_prefix = output_dir * "/neverworld_sixteenth"

arch   = GPU()
old_degree = 1
new_degree = 1/16

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

grid = MultiRegionGrid(grid.underlying_grid, partition = XPartition(2), devices = 2)

# Extend the vertical advection scheme
interp_init = false
init_file   = nothing #"files_lowres_new_bathy/restart_file_15_years.jld2" 

# Simulation parameters
Δt        = 0.5minutes
stop_time = 7000days

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

increase_simulation_Δt!(simulation, cutoff_time = 500days,  new_Δt = 2.0minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

checkpoint_time = 20days
checkpoint_outputs!(simulation, output_prefix; checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
