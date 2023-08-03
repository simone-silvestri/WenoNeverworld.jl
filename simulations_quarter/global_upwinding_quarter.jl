using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using JLD2

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/storage2/"
@show output_prefix = output_dir * "WenoNeverworldData/weno_fourth"

arch = GPU()
new_degree = 1/4
old_degree = 1/4

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70), H=7)
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70), H=7)
# orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70)) for old_degree = 1/4

# Extend the vertical advection scheme
interp_init = false
#init_file = "/home/sandre/Repositories/WenoNeverworld.jl/simulations_quarter/weno_four_checkpoint.jld2"
init_file = "/storage2/WenoNeverworldData/weno_fourth.jld2" #weno_four_checkpoint_iteration995943.jld2"

# Simulation parameters
Δt        = 20minutes
final_Δt  = 20minutes 
stop_time = 2000years

tracer_advection      = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(order = 9), 
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
coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme())
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          biharmonic_viscosity, free_surface, orig_grid, coriolis)

# Increase the time step up to final_Δt by final_day with "divisions" steps
starting_day = 0days
final_day = days # 365days
divisions = 5 # can only do it up to 5 times
Δts = range(Δt, final_Δt, length = divisions)
cutoff_times = range(starting_day, final_day, length = divisions)                                       
for (dt, cutoff_time) in zip(Δts, cutoff_times)
    println("Increasing Δt to $(prettytime(dt)) at $(prettytime(cutoff_time)) ")
    if dt > Δt
        increase_simulation_Δt!(simulation, cutoff_time = cutoff_time, new_Δt  = dt)
    end
end
#=
increase_simulation_Δt!(simulation, cutoff_time = 6 * 60days + 30days, new_Δt  = 2minutes)
increase_simulation_Δt!(simulation, cutoff_time = 6 * 60days + 60days, new_Δt  = 3minutes)
increase_simulation_Δt!(simulation, cutoff_time = 6 * 60days + 90days, new_Δt  = 4minutes)
increase_simulation_Δt!(simulation, cutoff_time = 6 * 60days + 120days, new_Δt = 5minutes)
=#
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
# standard_outputs!(simulation, output_prefix; overwrite_existing)
checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = years)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)

