using WenoNeverworld
using WenoNeverworld.Parameterizations
using WenoNeverworld.Auxiliaries
using Oceananigans
using Oceananigans.Units

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = GPU()

# The resolution in degrees
degree_resolution = 1/4

grid = NeverworldGrid(degree_resolution; arch)

# Do we need to interpolate? (interp_init) If `true` from which file?
interp_init = false # If interpolating from a different grid: `interp_init = true`
init_file   = nothing # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 10minutes
stop_time = 200years

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

# Here we use test `NNbackscatteringClosure` closure
horizontal_closure = NNbackscatteringClosure(eltype(grid); 
                                             architecture = arch, 
                                             weight_path = "model_weights.jld2")

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              wind_stress,
                                              buoyancy_relaxation,
                                              interp_init,
                                              horizontal_closure,
                                              init_file)
                                              
# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix)

# Initialize with a small time step and increase it after the 
# inital spin up has completed
increase_simulation_Δt!(simulation; cutoff_time = 20days,  new_Δt =  2minutes)
increase_simulation_Δt!(simulation; cutoff_time = 50days,  new_Δt =  4minutes)
increase_simulation_Δt!(simulation; cutoff_time = 100days, new_Δt =  6minutes)
increase_simulation_Δt!(simulation; cutoff_time = 150days, new_Δt =  8minutes)
increase_simulation_Δt!(simulation; cutoff_time = 200days, new_Δt = 10minutes)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

