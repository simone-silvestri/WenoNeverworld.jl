using WenoNeverworld
using WenoNeverworld: initial_buoyancy_parabola
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: HorizontalScalarBiharmonicDiffusivity

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_high_resolution"

arch   = GPU()
new_degree = 1/4

using CUDA
CUDA.device!(1)

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))

interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 5minutes
stop_time = 20years

using WenoNeverworld: geometric_νhb

momentum_advection   = VectorInvariant()
biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days) 

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, initial_buoyancy = initial_buoyancy_parabola, λ_buoy = 7days, momentum_advection, biharmonic_viscosity)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 1year
snapshot_time   = 1year
surface_time    = 30days
average_time    = 1year
standard_outputs!(simulation, output_prefix; checkpoint_time, snapshot_time, surface_time, average_time)

# initializing the time for wall_time calculation
run_simulation!(simulation)
