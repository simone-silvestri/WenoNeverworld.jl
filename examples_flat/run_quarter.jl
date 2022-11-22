using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity

output_dir    = joinpath(@__DIR__, "../files_four_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

interp_init = true
init_file   = "files_lowres_new_bathy/neverworld_lowres_checkpoint_iteration2067840.jld2" 

# Simulation parameters
Δt        = 0.5minutes
stop_time = 7000days

biharmonic_viscosity = HorizontalDivergenceScalarDiffusivity(ν = 50.0)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, biharmonic_viscosity, vertical_diffusivity)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 100days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
