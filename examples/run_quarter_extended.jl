using WenoNeverworld
using WenoNeverworld: initial_buoyancy_parabola
using Oceananigans

output_dir    = joinpath(@__DIR__, "../files_four_extended")
@show output_prefix = output_dir * "/neverworld_quarter_exetended"

arch   = GPU()
new_degree = 1/4

grid = NeverworldGridExtended(arch, new_degree)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = nothing

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

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
run_simulation!(simulation; init, init_file)
