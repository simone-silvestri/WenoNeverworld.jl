using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge, bathymetry_without_ridge

output_dir    = joinpath(@__DIR__, "../files_lowres_new_bathy")
@show output_prefix = output_dir * "/neverworld_lowres"

arch       = GPU()
new_degree = 1

orig_grid = NeverworldGrid(arch, new_degree; longitude = (-5, 65), bathymetry = bathymetry_with_ridge)
grid      = NeverworldGrid(arch, new_degree; longitude = (-5, 65))

# Remember to pass init file if we want to interpolate!
interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration2341318.jld2"

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 20minutes
stop_time = 120years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file) 

# Add outputs
checkpoint_time = 1year
snapshot_time   = 1year
surface_time    = 30days
average_time    = 1year
standard_outputs!(simulation, output_prefix; checkpoint_time, snapshot_time, surface_time, average_time, overwrite_existing=true)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
