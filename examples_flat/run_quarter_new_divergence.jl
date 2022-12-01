using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb

output_dir    = joinpath(@__DIR__, "../files_four_divergence_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter_new_divergence"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = false
init_file   = "files_four_new_bathy/neverworld_quarter_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 10minutes
stop_time = 100years

include("../src/vector_invariant_divergence.jl")

tracer_advection   = WENO(grid.underlying_grid) 
momentum_advection = WENO(vector_invariant = DivergenceStencil())

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, tracer_advection, momentum_advection)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
