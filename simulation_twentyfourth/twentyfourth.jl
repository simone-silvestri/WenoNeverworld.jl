using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using JLD2

# output_dir    = joinpath(@__DIR__, "./")
# @show output_prefix = output_dir * "/weno_two"
# output_dir = "/storage2/"
output_dir = "/home/sandre/Repositories/WenoNeverworld.jl/" # "/nobackup1/users/sandre/WenoNeverworldData/"
output_dir = "/pool001/users/sandre/WenoNeverworldData/"
@show output_prefix = output_dir * "weno_twentyfourth"

arch = GPU()
new_degree = 1/24
old_degree = 1/24

grid = NeverworldGrid(arch, new_degree, latitude = (-70, -20))
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, -20))

# Extend the vertical advection scheme
interp_init = false
# init_file = "/storage2/WenoNeverworldData/weno_eight_checkpoint_iteration648000.jld2"
# init_file = "/storage2/WenoNeverworldData/weno_eight_checkpoint_iteration142782.jld2"
init_file = "/pool001/users/sandre/WenoNeverworldData/weno_twentyfourth_checkpoint_iteration529462.jld2"

# Simulation parameters
Δt       =  10
final_Δt =  10 
stop_time = 200years

tracer_advection      = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(),
                                       divergence_scheme = WENO(),
                                         vertical_scheme = WENO(grid.underlying_grid))

biharmonic_viscosity  = nothing
convective_adjustment  = RiBasedVerticalDiffusivity()

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)
  
   return Int(ceil(2 * Δt / (CFL / wave_speed * local_Δ)))
end
  
free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(final_Δt, grid, g_Earth))
  
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          biharmonic_viscosity, free_surface, orig_grid)

# Increase the time step up to final_Δt by final_day with "divisions" steps
starting_day = 0days
final_day = 365days # 365days
divisions = 10 # can only do it up to 10 times
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

