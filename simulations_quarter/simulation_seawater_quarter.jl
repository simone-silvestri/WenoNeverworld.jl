using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using JLD2

output_prefix = "/nobackup/users/sbishnu/neverworld_seawater_quarter"

arch = GPU()
degree = 1/4

grid = NeverworldGrid(arch, degree, latitude = (-70, 70), H = 7)

interp_init = false
init_file   = nothing 

# Simulation parameters
max_Δt = 20minutes
min_Δt = 1minute
stop_time = 400years

tracer_advection      = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(order = 9),
                                         vertical_scheme = WENO(grid.underlying_grid))

free_surface = SplitExplicitFreeSurface(; grid, cfl = 0.75)
# Construct the neverworld simulation
simulation = neverworld_simulation_seawater(; grid, Δt = min_Δt, stop_time, interp_init, init_file,
                                              tracer_advection, momentum_advection,
                                              free_surface)

wizard = TimeStepWizard(; cfl = 0.3, max_Δt, min_Δt, max_change = 1.1)
simlation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
checkpoint_time = 1years
average_time    = 1years
snapshot_time   = 10years
checkpoint_outputs!(simulation, output_prefix; overwrite_existing, checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
