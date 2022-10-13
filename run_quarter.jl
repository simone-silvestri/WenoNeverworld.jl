using CUDA: device!
using CUDA
device!(3)

include("neverworld_utils.jl")

output_prefix = "neverworld_quarter"

H = 3

arch   = GPU()
old_degree = 1/4
new_degree = 1/4
grid      = neverworld_grid(arch, new_degree; H)
orig_grid = neverworld_grid(arch, old_degree; H)

init      = false

Δt        = 6minutes
stop_time = 100years

b_init = jldopen("evolved_initial_conditions_quarter.jld2")["b"]
u_init = jldopen("evolved_initial_conditions_quarter.jld2")["u"]
v_init = jldopen("evolved_initial_conditions_quarter.jld2")["v"]
w_init = jldopen("evolved_initial_conditions_quarter.jld2")["w"]
η_init = jldopen("evolved_initial_conditions_quarter.jld2")["η"]

if !(grid == orig_grid)
    @info "interpolating b field"
    b_init = interpolate_per_level(b_init, old_degree, new_degree, (Center, Center, Center), H)
    @info "interpolating u field"
    u_init = interpolate_per_level(u_init, old_degree, new_degree, (Face, Center, Center), H)
    @info "interpolating v field"
    v_init = interpolate_per_level(v_init, old_degree, new_degree, (Center, Face, Center), H)
    @info "interpolating w field"
    w_init = interpolate_per_level(w_init, old_degree, new_degree, (Center, Center, Face), H)
end

include("weno_neverworld.jl")
