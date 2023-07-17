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
old_degree = 1/4
new_degree = 1/2

orig_grid = NeverworldGrid(arch, old_degree) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = true
init_file   = "files_four_new_bathy/neverworld_quarter_checkpoint_iteration2763401.jld2" 

# Simulation parameters
Δt        = 1minutes
stop_time = 7000days

tracer_advection = WENO(grid.underlying_grid)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, tracer_advection)

increase_simulation_Δt!(simulation, cutoff_time = 300days,  new_Δt = 2minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
