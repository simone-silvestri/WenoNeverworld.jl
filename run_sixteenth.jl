include("initial_conditions.jl")
include("neverworld_utils.jl")

output_prefix = "files_sixteen/neverworld_sixteenth"

H = 5

arch   = GPU()
old_degree = 1/16
new_degree = 1/16

grid      = neverworld_grid(arch, new_degree; H)
orig_grid = neverworld_grid(arch, old_degree; H)

# initialize from scratch (or interpolated) - true, from file - false
init = true

# interpolate from old coarser solution - true (in combination with init = true)
interp_init = false

# file to initialize the simulation with or interpolate 
init_file = "files_four/neverworld_quarter_checkpoint_iteration96120.jld2" 

Δt        = 1minutes
stop_time = 0.5years

if interp_init
    b_init = jldopen(init_file)["b"][H+1:end-H, H+1:end-H, H+1:end-H]
    u_init = jldopen(init_file)["u"][H+1:end-H, H+1:end-H, H+1:end-H]
    v_init = jldopen(init_file)["v"][H+1:end-H, H+1:end-H, H+1:end-H]
    η_init = jldopen(init_file)["η"][H+1:end-H, H+1:end-H, H+1:end-H]
    if !(grid == orig_grid)
         @info "interpolating b field"
         b_init = interpolate_per_level(b_init, old_degree, new_degree, (Center, Center, Center), H)
         @info "interpolating u field"
         u_init = interpolate_per_level(u_init, old_degree, new_degree, (Face, Center, Center), H)
         @info "interpolating v field"
         v_init = interpolate_per_level(v_init, old_degree, new_degree, (Center, Face, Center), H)
         @info "interpolating w field"
         η_init = interpolate_per_level(η_init, old_degree, new_degree, (Center, Center, Face), H)
    end
end

include("weno_neverworld.jl")

function increase_Δt(simulation)
    if simulation.model.clock.time > 90days
        simulation.Δt = 2minutes
    end
end

simulation.callbacks[:change_Δt] = Callback(increase_Δt, TimeInterval(90days))

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

if init
    run!(simulation)
else
    clock = jldopen(init_file)["clock"]
    model.clock.time      = clock.time
    model.clock.iteration = clock.iteartion

    run!(simulation, pickup=init_file)
end

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
