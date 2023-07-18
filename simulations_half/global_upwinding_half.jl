using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using JLD2

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/weno_quarter"

arch = GPU()
new_degree = 1/2

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))

# Extend the vertical advection scheme
interp_init = false
init_file   = nothing # "./global_upwinding_quarter_weno_checkpoint_iteration3020910.jld2"

# Simulation parameters
Δt        = 3minutes
stop_time = 200years

tracer_advection      = WENO(grid.underlying_grid) 
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(), 
                                       divergence_scheme = WENO(),
                                         vertical_scheme = WENO(grid.underlying_grid))

biharmonic_viscosity  = nothing

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / minimum_xspacing(grid)^2 + 1 / minimum_yspacing(grid)^2)

    return Int(ceil(2 * Δt / (CFL / wave_speed * local_Δ)))
end

free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(10minutes, grid, g_Earth))

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          free_surface)

initial_buoyancy_one(x, y, z) = WenoNeverworld.exponential_profile(z)

set!(simulation.model.tracers.b, initial_buoyancy_one)

increase_simulation_Δt!(simulation, cutoff_time =  30days, new_Δt =  5minutes)
increase_simulation_Δt!(simulation, cutoff_time =  60days, new_Δt =  8minutes)
increase_simulation_Δt!(simulation, cutoff_time =  120days, new_Δt = 10minutes)
increase_simulation_Δt!(simulation, cutoff_time =  180days, new_Δt = 15minutes)
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
checkpoint_time = 1years
average_time    = 1years
snapshot_time   = 10years
checkpoint_outputs!(simulation, output_prefix; overwrite_existing, checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
