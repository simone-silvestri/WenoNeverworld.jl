using Oceananigans
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../files_four_centered_divergence")
@show output_prefix = output_dir * "/neverworld_quarter_centered_divergence"

arch   = GPU()
old_degree = 1/4
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree)
grid      = NeverworldGrid(arch, new_degree)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = "files_four_centered/neverworld_quarter_centered_checkpoint_iteration1066000.jld2"

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 6minutes
stop_time = 20.5years

## Changing parameterizations
using WenoNeverworld: geometric_νhb

biharmonic_viscosity = HorizontalDivergenceScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)
momentum_advection   = VectorInvariant()

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, biharmonic_viscosity, momentum_advection, Δt, stop_time, interp_init, init_file)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 30days
standard_outputs!(simulation, output_prefix; checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
