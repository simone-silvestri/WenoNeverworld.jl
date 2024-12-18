using WenoNeverworld
using WenoNeverworld.NeverworldGrids
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using CUDA

CUDA.device!(3)
arch = GPU()

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/storage4/WenoNeverworldData/eighth_degree_new/"
@show output_prefix = output_dir * "weno_eighth"

# The resolution in degrees
new_degree = 1/8
old_degree = 1/4

z_faces = exponential_z_faces(; Nz = 35, depth = 3000)
grid = NeverworldGrid(new_degree; arch, z_faces)
previous_grid = NeverworldGrid(old_degree; arch, z_faces)


# Do we need to interpolate? (interp_init) If `true` from which file?
interp_init = true # If interpolating from a different grid: `interp_init = true`
init_file   = "/storage4/WenoNeverworldData/quarter_degree_new/weno_quarter__checkpoint_iteration70005600.jld2" # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 5minutes
stop_time = 4000years

# Latitudinal wind stress acting on the zonal velocity
# a piecewise-cubic profile interpolated between
# x = φs (latitude) and y = τs (stress)
φs = (-70.0, -45.0, -15.0,  0.0,  15.0, 45.0, 70.0)
τs = (  0.0,   0.2,  -0.1, -0.02, -0.1,  0.1,  0.0)
wind_stress = WindStressBoundaryCondition(; φs, τs)

# Buoyancy relaxation profile:
# a parabolic profile between 0, at the poles, and ΔB = 0.06 at the equator
# the restoring time is λ = 7days
buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(ΔB = 0.06, λ = 7days)

# Wanna use a different profile? Try this:
# @inline seasonal_cosine_scaling(y, t) = cos(π * y / 70) * sin(2π * t / 1year)
# buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(seasonal_cosine_scaling; ΔB = 0.06, λ = 7days)    

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; previous_grid, Δt, stop_time,
                                              wind_stress,
                                              buoyancy_relaxation,
                                              interp_init,
                                              init_file, vertical_diffusivity = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=3e-5))
                                              
model = simulation.model

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 30days)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

