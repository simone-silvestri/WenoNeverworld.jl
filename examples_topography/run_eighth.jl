using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using WenoNeverworld: bathymetry_with_ridge

output_dir    = joinpath(@__DIR__, "../files_eight")
@show output_prefix = output_dir * "/neverworld_eighth"

arch   = GPU()
old_degree = 1
new_degree = 1/8

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65), bathymetry = bathymetry_with_ridge)
grid      = NeverworldGrid(arch, new_degree; bathymetry = bathymetry_with_ridge)

interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration3313678.jld2"

# Simulation parameters
Δt        = 2minutes
stop_time = 20years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

increase_simulation_Δt!(simulation, cutoff_time = 30days,  new_Δt = 2.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 60days,  new_Δt = 3minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
