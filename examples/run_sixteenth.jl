using Oceananigans
using Oceananigans.Units
using WenoNeverworld

output_prefix = joinpath(@__DIR__, "/neverworld_sixteenth")

H = 5

arch   = GPU()
old_degree = 1/4
new_degree = 1/16

grid      = neverworld_grid(arch, new_degree; H)
orig_grid = neverworld_grid(arch, old_degree; H)

# Remember to pass init file if we want to interpolate!
interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 1minutes
stop_time = 0.5years

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file)

# Increase simulation Δt after 50days
increase_simulation_Δt!(simulation, cutoff_time = 50days, new_Δt = 2minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

# Add outputs
checkpoint_time = 50days
standard_outputs!(simulation, output_prefix; checkpoint_time)

if init
    run!(simulation)
else
    update_simulation_clock!(simulation, init_file)
    run!(simulation, pickup=init_file)
end

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
