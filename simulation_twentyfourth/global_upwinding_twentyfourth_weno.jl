using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.Distributed
using JLD2

using MPI
MPI.Init()

rank   = MPI.Comm_rank(MPI.COMM_WORLD)
Nranks = MPI.Comm_size(MPI.COMM_WORLD)

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/weno_sixteenth"

child_arch = GPU()

arch = DistributedArch(child_arch, ranks = (Nranks, 1, 1), use_buffers = true)

new_degree = 1/24
old_degree = 1/4

grid      = NeverworldGrid(arch, new_degree, latitude = (-70, -20))
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, -20))

# Extend the vertical advection scheme
interp_init = true
init_file   = "/home/sandre/Repositories/WenoNeverworld.jl/simulations_quarter/weno_four_checkpoint_iteration46346.jld2"

# Simulation parameters
Δt        = 20
stop_time = 100years

final_Δt  = 90 

tracer_advection      = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(),
                                       divergence_scheme = WENO(),
                                         vertical_scheme = WENO(grid.underlying_grid))

biharmonic_viscosity  = nothing
convective_adjustment  = RiBasedVerticalDiffusivity()

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / minimum_xspacing(grid)^2 + 1 / minimum_yspacing(grid)^2)
  
   return Int(ceil(2 * Δt / (CFL / wave_speed * local_Δ)))
end
  
free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(final_Δt, grid, g_Earth))
  
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          biharmonic_viscosity, free_surface, orig_grid)

starting_day = 0days
final_day    = 100days # 365days
divisions    = 5 # can only do it up to 5 times
Δts = range(Δt, final_Δt, length = divisions)
cutoff_times = range(starting_day, final_day, length = divisions)                                       

for (dt, cutoff_time) in zip(Δts, cutoff_times)
    println("Increasing Δt to $(prettytime(dt)) at $(prettytime(cutoff_time)) ")
    if dt > Δt
        increase_simulation_Δt!(simulation, cutoff_time = cutoff_time, new_Δt  = dt)
    end
end

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)

