using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb

output_dir    = joinpath(@__DIR__, "../files_eight_new_bathy")
@show output_prefix = output_dir * "/neverworld_eighth"

arch   = GPU()
old_degree = 1
new_degree = 1/8

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = true
init_file   = "files_lowres_new_bathy/neverworld_lowres_checkpoint_iteration2067840.jld2" 

# Simulation parameters
Δt        = 2minutes
stop_time = 7000days

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

increase_simulation_Δt!(simulation, cutoff_time = 300days,  new_Δt = 3.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 600days,  new_Δt = 5.0minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
