using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using WenoNeverworld: geometric_νhb
using Oceananigans.BuoyancyModels: g_Earth

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/global_upwinding_sixteen_leith"

arch = GPU()
new_degree = 1/16

grid = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = true
init_file   = "/pool001/users/ssilvest/old_weno_sixteen/global_upwinding_sixteen_leith_checkpoint_iteration1861600.jld2"

# Simulation parameters
Δt        = 2minutes
stop_time = 100years

using WenoNeverworld: leith_laplacian_viscosity

flux_form_weno        = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant()

biharmonic_viscosity  = leith_laplacian_viscosity(C_vort = 2.0, C_div = 2.0)

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / minimum_xspacing(grid)^2 + 1 / minimum_yspacing(grid)^2)
  
   return Int(ceil(2 * Δt / (CFL / wave_speed * local_Δ)))
end
  
free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(Δt, grid, g_Earth))
  
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection = flux_form_weno, momentum_advection, 
                                          biharmonic_viscosity, free_surface)

increase_simulation_Δt!(simulation, cutoff_time = 30days, new_Δt = 2minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
checkpoint_time = 1years
standard_outputs!(simulation, output_prefix; checkpoint_time, overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)

