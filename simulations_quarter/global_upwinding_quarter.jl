using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using JLD2

years = 43200 * 365 # every half a year
output_dir    = joinpath(@__DIR__, "./")
output_dir = "/pool001/users/sandre/WenoNeverworldData/"
output_dir = "/orcd/nese/raffaele/001/sandre/WenoNeverworld/"
@show output_prefix = output_dir * "weno_fourth"

arch = GPU()
new_degree = 1/4
old_degree = 1/4

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70), H = 7)
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70), H = 7)
# orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70)) for old_degree = 1/4

# Extend the vertical advection scheme
interp_init = false
# init_file = "/pool001/users/sandre/WenoNeverworldData/weno_fourth_checkpoint_iteration19606397.jld2"
# init_file = "/storage2/WenoNeverworldData/weno_four_checkpoint_iteration2630343.jld2"
init_file = "/orcd/nese/raffaele/001/sandre/WenoNeverworld/weno_fourth.jld2"

# Simulation parameters, can probably do 16 minutes, previously was doing 12.5 minutes
Δt       = 16minutes
final_Δt = 16minutes 
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
final_day = 365days # 365days
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

#=
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
checkpoint_time = 1years
average_time    = 1years
snapshot_time   = 10years
checkpoint_outputs!(simulation, output_prefix; overwrite_existing, checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
=#

