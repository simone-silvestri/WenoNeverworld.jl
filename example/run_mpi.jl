using MPI
MPI.Init()

using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.Distributed

output_dir    = joinpath(@__DIR__, "./")
output_dir = "./"
@show output_prefix = output_dir * "WenoNeverworldData/weno_two" 

rx = parse(Int, get(ENV, "RX", "1"))
ry = parse(Int, get(ENV, "RY", "1"))

arch = DistributedArch(GPU(), ranks = (rx, ry, 1), topology = (Periodic, Bounded, Bounded))

# The resolution in degrees
degree = 1 / 48 # 1 / 64 degree resolution

grid = NeverworldGrid(degree; arch)

# Extend the vertical advection scheme
interp_init = false # Do we need to interpolate? (interp_init) If `true` from which file? # If interpolating from a different grid: `interp_init = true`
init_file   = nothing # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 1minutes
stop_time = 3000years

free_surface = SplitExplicitFreeSurface(; grid, cfl = 0.75, fixed_Δt = 2minutes)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time, 
                                              free_surface,
                                              interp_init,
                                              init_file)
                                              

# Adaptable time step
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 45minutes, min_Δt = 2minutes, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
# checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)
# vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)
