using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization

output_dir    = joinpath(@__DIR__, "../files_sixteen")
@show output_prefix = output_dir * "/neverworld_sixteenth"

arch   = GPU()
old_degree = 1
new_degree = 1/16

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65))
grid      = NeverworldGrid(arch, new_degree)

interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration3313678.jld2"

# Simulation parameters
Δt        = 1minutes
stop_time = 2years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

increase_simulation_Δt!(simulation, cutoff_time = 30days,  new_Δt = 2minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
