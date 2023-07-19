using WenoNeverworld
using WenoNeverworld: initial_buoyancy_parabola
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = GPU()
degree_resolution = 1/4

grid = NeverworldGrid(arch, degree_resolution)

interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 10minutes
stop_time = 200years

# Construct the neverworld simulation
simulation = neverworld_simulation_seawater(; grid, Δt, stop_time) 



# Add outputs
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

