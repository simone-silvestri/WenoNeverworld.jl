using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/uq_one_degree_temperature"

arch = GPU()
new_degree = 1

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))

# Pickup data
interp_init = true
pickup_data = true
init_file   = "/nobackup/users/hannahlu/WenoNeverworld.jl_old/examples_full_latitude/outputs_data/neverworld_high_resolution_checkpoint_iteration5256478.jld2"

# Simulation parameters
Δt        = 20minutes
stop_time = 30years

using WenoNeverworld: geometric_νhb

tracer_advection     = WENO(grid.underlying_grid) 
momentum_advection   = VectorInvariant()
biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 2days)

# Parameters
τₚ = [0.2, -0.1, -0.02, -0.1, 0.1]
Sₚ = [-2e-8, 2e-8, -4e-8, 2e-8, -2e-8]
ΔT = 35.0

# Construct the neverworld simulation
simulation = neverworld_simulation_seawater(; grid, Δt, stop_time, interp_init, init_file, 
                                              tracer_advection, momentum_advection, 
                                              biharmonic_viscosity,
                                              τₚ, Sₚ, ΔT, pickup_data)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
checkpoint_time = 1years
average_time    = 1years
snapshot_time   = 10years
standard_seawater_outputs!(simulation, output_prefix; overwrite_existing, checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
