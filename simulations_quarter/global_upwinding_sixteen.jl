using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
#using CairoMakie # You have to add this to your global enviroment: `] add CairoMakie`

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/WenoNeverworldData"
@show output_prefix = output_dir * "weno_sixteen"

arch = GPU()

# The resolution in degrees
degree_resolution = 1/16
new_degree = 1/16
old_degree = 1/16

grid = NeverworldGrid(new_degree; arch)
previous_grid = NeverworldGrid(old_degree; arch)

# Extend the vertical advection scheme
interp_init = true # Do we need to interpolate? (interp_init) If `true` from which file? # If interpolating from a different grid: `interp_init = true`
init_file = "/WenoNeverworldData/weno_sixteenth_checkpoint_iteration1826520.jld2" # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 10minutes
stop_time = 3000years

# Latitudinal wind stress acting on the zonal velocity
# a piecewise-cubic profile interpolated between
# x = φs (latitude) and y = τs (stress)
#φs = (-70.0, -45.0, -15.0,  0.0,  15.0, 45.0, 70.0)
#τs = (  0.0,   0.2,  -0.1, -0.02, -0.1,  0.1,  0.0)
#wind_stress = WindStressBoundaryCondition(; φs, τs)

# Buoyancy relaxation profile:
# a parabolic profile between 0, at the poles, and ΔB = 0.06 at the equator
# the restoring time is λ = 7days
#buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(ΔB = 0.06, λ = 7days)

# Wanna use a different profile? Try this:
# @inline seasonal_cosine_scaling(y, t) = cos(π * y / 70) * sin(2π * t / 1year)
# buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(seasonal_cosine_scaling; ΔB = 0.06, λ = 7days)    

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time, 
                                              interp_init,
                                              init_file)
                                              

# Adaptable time step
# wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 40minutes, min_Δt = 15minutes, max_change = 1.1)
# simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)
#vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)