using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb

output_dir    = joinpath(@__DIR__, "../files_four_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = true
init_file   = "files_lowres_new_bathy/restart_file_15_years.jld2" 

# Simulation parameters
Δt        = 2minutes
stop_time = 7000days

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 5.0minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 300days, new_Δt = 10minutes)
increase_simulation_Δt!(simulation, cutoff_time = 400days, new_Δt = 15minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
