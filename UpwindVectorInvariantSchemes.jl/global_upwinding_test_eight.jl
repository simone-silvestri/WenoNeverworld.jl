using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity, HorizontalDivergenceFormulation
using WenoNeverworld: geometric_νhb
using Oceananigans.Advection: VelocityStencil
using UpwindVectorInvariantSchemes

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/global_upwinding_eight"

using CUDA
CUDA.device!(3)

arch = GPU()
old_degree = 1/4
new_degree = 1/8

orig_grid = NeverworldGrid(arch, old_degree)
grid     = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = false
init_file   = "./global_upwinding_eight_checkpoint_iteration382000.jld2" 

# Simulation parameters
Δt        = 1minutes
stop_time = 100years

flux_form_weno        = WENO(grid.underlying_grid) 
vector_invariant_weno = WENO(VelocityStencil())
momentum_advection    = GlobalVectorInvariant(; ζ_upwind_scheme = vector_invariant_weno, vertical_scheme = WENO(grid.underlying_grid))

convective_adjustment  = RiBasedVerticalDiffusivity()
biharmonic_viscosity   = nothing

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, tracer_advection = flux_form_weno, momentum_advection, biharmonic_viscosity, convective_adjustment)

increase_simulation_Δt!(simulation, cutoff_time = 30days, new_Δt = 2minutes)
# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = false
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
