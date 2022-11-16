using Oceananigans
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../files_sixteen")
@show output_prefix = output_dir * "/neverworld_sixteenth"

H = 5

arch   = GPU()
old_degree = 1/4
new_degree = 1/16

orig_grid = neverworld_grid(arch, old_degree; H)
grid      = neverworld_grid(arch, new_degree; H)

# Remember to pass init file if we want to interpolate!
interp_init = true
init_file   = "files_eight/neverworld_eighth_checkpoint_iteration422900.jld2"

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 1minutes
stop_time = 0.5years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

# Increase simulation Δt after 50days
increase_simulation_Δt!(simulation, cutoff_time = 0.3years, new_Δt = 2minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 20days
standard_outputs!(simulation, output_prefix; checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)