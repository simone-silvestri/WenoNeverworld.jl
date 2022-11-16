using Oceananigans
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../files_eight")
@show output_prefix = output_dir * "/neverworld_eighth"

arch   = GPU()
old_degree = 1/4
new_degree = 1/8

orig_grid = NeverworldGrid(arch, old_degree)
grid      = NeverworldGrid(arch, new_degree)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = "files_eight/neverworld_eighth_checkpoint_iteration422900.jld2"

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 3minutes
stop_time = 5years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

# # Increase simulation Δt after 40days
# increase_simulation_Δt!(simulation, cutoff_time = 150days, new_Δt = 3minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 0.5years
standard_outputs!(simulation, output_prefix; checkpoint_time, overwrite_existing = false)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
