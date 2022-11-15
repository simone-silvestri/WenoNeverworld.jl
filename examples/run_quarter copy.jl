using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../files_four_centered")
@show output_prefix = output_dir * "/neverworld_quarter"

H = 5

arch   = GPU()
old_degree = 1/4
new_degree = 1/4

orig_grid = neverworld_grid(arch, old_degree; H)
grid      = neverworld_grid(arch, new_degree; H)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = nothing
init        = true

# Simulation parameters
Δt        = 6minutes
stop_time = 20years

using WenoNeverworld: geometric_νhb

biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(geometric_νhb, discrete_form=true, parameters = 5days)
momentum_advection   = VectorInvariant()

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, biharmonic_viscosity, momentum_advection, Δt, stop_time, interp_init, init_file)

# Increase simulation Δt after 40days
increase_simulation_Δt!(simulation, cutoff_time = 150days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 1year
standard_outputs!(simulation, output_prefix; checkpoint_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
