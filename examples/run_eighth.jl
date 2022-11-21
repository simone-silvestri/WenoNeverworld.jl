using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization

output_dir    = joinpath(@__DIR__, "../files_eight")
@show output_prefix = output_dir * "/neverworld_eighth"

arch   = GPU()
old_degree = 1
new_degree = 1/8

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65))
grid      = NeverworldGrid(arch, new_degree)

# Remember to pass init file if we want to interpolate!
interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration3313678.jld2"

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 2minutes
stop_time = 20years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, vertical_diffusivity, convective_adjustment)

increase_simulation_Δt!(simulation, cutoff_time = 30days,  new_Δt = 2.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 60days,  new_Δt = 3minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
