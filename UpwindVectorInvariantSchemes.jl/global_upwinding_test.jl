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
@show output_prefix = output_dir * "/global_upwinding"

arch = GPU()
new_degree = 1/4

grid = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = false
init_file   = "../files_four_new_bathy/neverworld_quarter_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 10minutes
stop_time = 100years

tracer_advection = WENO(grid.underlying_grid) 
upwind_scheme    = WENO(VelocityStencil())
momentum_advection = GlobalVectorInvariant(; upwind_scheme, vertical_scheme = tracer_advection)

biharmonic_viscosity = nothing

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, tracer_advection, momentum_advection, biharmonic_viscosity)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
