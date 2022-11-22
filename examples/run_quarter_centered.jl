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

grid      = NeverworldGrid(arch, new_degree; bathymetry = bathymetry_with_ridge)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = "files_four_centered/restart_file.jld2"

using WenoNeverworld: geometric_νhb

biharmonic_viscosity  = HorizontalScalarBiharmonicDiffusivity(ν=geometric_νhb, discrete_form=true, parameters = 5days)
momentum_advection    = VectorInvariant()

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 5minutes
stop_time = 40years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, biharmonic_viscosity, momentum_advection, μ_drag = 0.003)

increase_simulation_Δt!(simulation, cutoff_time = 7300days,  new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 7400days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
