using MPI
MPI.Init()

using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture, minimum_xspacing, minimum_yspacing
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: all_reduce
using Oceananigans.Models.HydrostaticFreeSurfaceModels: FixedSubstepNumber

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/nobackup/users/lcbrock/WenoNeverworldData/"
@show output_prefix = output_dir * "weno_thirtytwo" 

Rx = parse(Int, get(ENV, "RX", "1"))
Ry = parse(Int, get(ENV, "RY", "1"))

arch = Distributed(GPU(), partition = Partition(Rx, Ry))

# The resolution in degrees
degree = 1 / 32 # degree resolution
previous_degree = 1 /8 

grid = NeverworldGrid(degree; arch)
# previous_grid needs to be on another architecture!!!
previous_grid = NeverworldGrid(previous_degree; arch = CPU())

# Extend the vertical advection scheme
interp_init = true # Do we need to interpolate? (interp_init) If `true` from which file? # If interpolating from a different grid: `interp_init = true`
init_file   = output_dir * "weno_eight_checkpoint_iteration15855224.jld2" # "test_mpi_" * string(MPI.Comm_rank(MPI.COMM_WORLD)) * "_checkpoint_iteration_" # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 0.01minutes
stop_time = 3000years
max_Δt    = 2minutes
using Oceananigans.Models.HydrostaticFreeSurfaceModels: FixedTimeStepSize
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing

substepping = FixedTimeStepSize(; cfl = 0.75, grid)
@show arch.local_rank, substepping.Δt_barotropic, grid.Lz, minimum_xspacing(grid), minimum_yspacing(grid)
substeps    = ceil(Int, 2*max_Δt / substepping.Δt_barotropic)
substeps    = all_reduce(max, substeps, arch)

free_surface = SplitExplicitFreeSurface(; substeps)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              previous_grid,
                                              free_surface,
                                              interp_init,
                                              init_file)
                                              

# Adaptable time step
wizard = TimeStepWizard(; cfl = 0.35, max_Δt, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 30days)
# vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)
