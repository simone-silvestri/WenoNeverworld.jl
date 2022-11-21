using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.Advection: VelocityStencil

output_dir    = joinpath(@__DIR__, "../files_four_highorder")
@show output_prefix = output_dir * "/neverworld_quarter_highorder"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65))
grid      = NeverworldGrid(arch, new_degree; H = 9)

# Remember to pass init file if we want to interpolate!
interp_init = true
init_file   = "files_lowres/neverworld_lowres_checkpoint_iteration3313678.jld2"

vertical_diffusivity = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν = 1e-4, κ = 1e-5)
convective_adjustment= ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.5)
momentum_advection   = WENO(order = 9, vector_invariant = VelocityStencil()) 

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 2minutes
stop_time = 20years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, vertical_diffusivity, convective_adjustment, momentum_advection)

increase_simulation_Δt!(simulation, cutoff_time = 60days,  new_Δt = 5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 90days,  new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 120days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
