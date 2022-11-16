using Oceananigans
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../files_two")
@show output_prefix = output_dir * "/neverworld_half"

arch   = GPU()
old_degree = 1/2
new_degree = 1/2

orig_grid = NeverworldGrid(arch, old_degree)
grid      = NeverworldGrid(arch, new_degree)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = nothing

# init always has to be true with interp_init, otherwise it depends if we start from a file or not
init = interp_init ? true : (init_file isa Nothing ? true : false)

# Simulation parameters
Δt        = 10minutes
stop_time = 100years

preconditioner_method   = :SparseInverse
preconditioner_settings = (ε = 0.001, nzrel = 5) 

free_surface         = ImplicitFreeSurface(; preconditioner_method, preconditioner_settings)
biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, free_surface, biharmonic_viscosity)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

increase_simulation_Δt(simualation, cutoff_time = 60days, new_Δt = 15minutes)
# Add outputs
checkpoint_time = 1year
snapshot_time   = 365days
surface_time    = 10days
average_time    = 365days
standard_outputs!(simulation, output_prefix; checkpoint_time, snapshot_time, surface_time, average_time)

# initializing the time for wall_time calculation
run_simulation!(simulation; init, init_file)
