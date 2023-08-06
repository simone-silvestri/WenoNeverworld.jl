using WenoNeverworld
using Oceananigans
using Oceananigans.Units

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = GPU()
degree_resolution = 1/4

grid = NeverworldGrid(degree_resolution; arch)

# Do we need to interpolate? (interp_init) If `true` from which file?
interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 10minutes
stop_time = 200years

# Construct the neverworld simulation
simulation = neverworld_simulation_seawater(grid; Δt, stop_time) 

# Maybe increase the time step to a new time-step
increase_simulation_Δt!(simulation, cutoff_time = 60days, new_Δt = 20minutes)

# Add outputs
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

