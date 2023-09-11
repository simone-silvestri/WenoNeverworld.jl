using Preferences
const iscray = parse(Bool, load_preference(Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"), "iscray", "false"))
@debug "Preloading GTL library" iscray
if iscray
    import Libdl
    Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
end

using MPI
MPI.Init()

using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.Distributed
#using CairoMakie # You have to add this to your global enviroment: `] add CairoMakie`

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/storage2/"
@show output_prefix = output_dir * "WenoNeverworldData/weno_two" 

rx = parse(Int, get(ENV, "RX", "1"))
ry = parse(Int, get(ENV, "RY", "1"))

arch = DistributedArch(GPU(), ranks = (rx, ry, 1), topology = (Periodic, Bounded, Bounded))

# The resolution in degrees
degree = 1 / 64 # 1 / 64 degree resolution

grid = NeverworldGrid(degree; arch)

# Extend the vertical advection scheme
interp_init = true # Do we need to interpolate? (interp_init) If `true` from which file? # If interpolating from a different grid: `interp_init = true`
init_file   = nothing # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 1minutes
stop_time = 3000years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; previous_grid, Δt, stop_time, 
                                              interp_init,
                                              init_file)
                                              

# Adaptable time step
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 45minutes, min_Δt = 15minutes, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
# checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)
# vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 10years)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)