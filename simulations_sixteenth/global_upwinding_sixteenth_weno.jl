using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using JLD2

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/weno_sixteenth"

arch = GPU()
new_degree = 1/16
old_degree = 1/4

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70))

# Extend the vertical advection scheme
interp_init = true
init_file   = "/home/sandre/Repositories/WenoNeverworld.jl/restart_files/weno_quarter_checkpoint.jld2" 

# Simulation parameters
Δt        = 1minutes
stop_time = 100years

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
  
free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(Δt, grid, g_Earth))
  
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          biharmonic_viscosity, free_surface, orig_grid)

increase_simulation_Δt!(simulation, cutoff_time = 30days, new_Δt = 2minutes)
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)

