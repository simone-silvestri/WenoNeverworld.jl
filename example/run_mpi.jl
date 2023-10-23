using MPI
MPI.Init()

using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: all_reduce

output_dir    = joinpath(@__DIR__, "./")
output_dir = "./"
@show output_prefix = output_dir * "test_mpi" 

Rx = parse(Int, get(ENV, "RX", "1"))
Ry = parse(Int, get(ENV, "RY", "1"))

arch = Distributed(GPU(), partition = Partition(Rx, Ry))

# The resolution in degrees
degree = 1 / 64 # degree resolution

grid = NeverworldGrid(degree; arch)

# Extend the vertical advection scheme
interp_init = false # Do we need to interpolate? (interp_init) If `true` from which file? # If interpolating from a different grid: `interp_init = true`
init_file   = nothing # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 0.01minutes
stop_time = 3000years
max_Δt    = 2minutes
using Oceananigans.Models.HydrostaticFreeSurfaceModels: FixedTimeStepSize

substepping = FixedTimeStepSize(; cfl = 0.75, grid)
substeps    = ceil(Int, 2*max_Δt / substepping.Δt_barotropic)
substeps    = all_reduce(max, substeps, arch)

free_surface = SplitExplicitFreeSurface(; substeps)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time, 
                                              free_surface,
                                              interp_init,
                                              init_file)
                                              

# Adaptable time step
wizard = TimeStepWizard(; cfl = 0.35, max_Δt, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
# checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)
# vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)
