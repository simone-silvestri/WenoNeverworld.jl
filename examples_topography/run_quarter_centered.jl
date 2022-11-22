using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization

output_dir    = joinpath(@__DIR__, "../files_four_centered")
@show output_prefix = output_dir * "/neverworld_quarter_centered"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65), bathymetry = bathymetry_with_ridge) 
grid      = NeverworldGrid(arch, new_degree; bathymetry = bathymetry_with_ridge)

# Remember to pass init file if we want to interpolate!
interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration3313678.jld2" 

using WenoNeverworld: cosine_νhb

biharmonic_viscosity  = HorizontalScalarBiharmonicDiffusivity(ν=cosine_νhb, discrete_form=true, parameters = 1e11)
momentum_advection    = VectorInvariant()

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 2minutes
stop_time = 7000days

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, biharmonic_viscosity, momentum_advection, μ_drag = 0.003)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 100days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
