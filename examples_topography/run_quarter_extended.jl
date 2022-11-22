using WenoNeverworld
using WenoNeverworld: initial_buoyancy_parabola
using Oceananigans
using Oceananigans.Units

output_dir    = joinpath(@__DIR__, "../files_four_extended")
@show output_prefix = output_dir * "/neverworld_quarter_exetended"

arch   = GPU()
new_degree = 1/4

grid = NeverworldGridExtended(arch, new_degree)

interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 5minutes
stop_time = 20years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, initial_buoyancy = initial_buoyancy_parabola)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 1year
snapshot_time   = 1year
surface_time    = 30days
average_time    = 1year
standard_outputs!(simulation, output_prefix; checkpoint_time, snapshot_time, surface_time, average_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
