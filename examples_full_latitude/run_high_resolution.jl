using WenoNeverworld
using WenoNeverworld: initial_buoyancy_parabola
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_high_resolution"

using CUDA
CUDA.device!(3)

arch        = GPU()
resolution = 1

grid = NeverworldGrid(arch, resolution, latitude = (-70, 70))

interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 20minutes
stop_time = 200years

using WenoNeverworld: geometric_νhb

momentum_advection   = VectorInvariant()
biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days) 

# Construct the neverworld simulation
simulation = neverworld_simulation_seawater(; grid, Δt, stop_time, momentum_advection, biharmonic_viscosity)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
