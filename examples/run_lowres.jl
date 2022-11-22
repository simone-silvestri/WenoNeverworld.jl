using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge, bathymetry_without_ridge

output_dir    = joinpath(@__DIR__, "../files_lowres_new_bathy")
@show output_prefix = output_dir * "/neverworld_lowres"

arch       = GPU()
new_degree = 1

grid = NeverworldGrid(arch, new_degree; longitude = (-5, 65))

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = nothing

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 2.5minutes
stop_time = 40years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file) 

increase_simulation_Δt!(simulation, cutoff_time = 30days, new_Δt = 5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 60days, new_Δt = 10minutes)

# Add outputs
checkpoint_time = 1year
snapshot_time   = 1year
surface_time    = 30days
average_time    = 1year
standard_outputs!(simulation, output_prefix; checkpoint_time, snapshot_time, surface_time, average_time, overwrite_existing=true)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
