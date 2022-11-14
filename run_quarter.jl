include("initial_conditions.jl")
include("neverworld_utils.jl")

output_prefix = "files_four/neverworld_quarter"

H = 5

arch   = GPU()
old_degree = 1/4
new_degree = 1/4

grid      = neverworld_grid(arch, new_degree; H)
orig_grid = neverworld_grid(arch, old_degree; H)

# initialize from scratch (or interpolated) - true, from file - false
init = false

# interpolate from old coarser solution - true (in combination with init = true)
interp_init = false

# file to initialize the simulation with or interpolate 
init_file = "files_four/neverworld_quarter_checkpoint_iteration2382480.jld2"

Δt        = 10minutes
stop_time = 45years
checkpoint_time = 0.5years

include("weno_neverworld.jl")

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

if init
    run!(simulation)
else
    clock = jldopen(init_file)["clock"]
    simulation.model.clock.time = clock.time
    simulation.model.clock.iteration = clock.iteration

    run!(simulation, pickup=init_file)
end

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
