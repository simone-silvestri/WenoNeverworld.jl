using CUDA: device!
device!(1)

include("neverworld_utils.jl")

H = 5
H_orig = 5
arch   = GPU()
old_degree = 1/16
new_degree = 1/16
grid      = neverworld_grid(arch, new_degree; H = H)
orig_grid = neverworld_grid(arch, old_degree; H = H_orig)

b_init = jldopen("evolved_initial_conditions_sixteenth.jld2")["b"]
u_init = jldopen("evolved_initial_conditions_sixteenth.jld2")["u"]
v_init = jldopen("evolved_initial_conditions_sixteenth.jld2")["v"]
w_init = jldopen("evolved_initial_conditions_sixteenth.jld2")["w"]
η_init = jldopen("evolved_initial_conditions_sixteenth.jld2")["η"]

if !(grid == orig_grid)
    @info "interpolating b field"
    b_init = interpolate_per_level(b_init, old_degree, new_degree, (Center, Center, Center), H, H_orig)
    @info "interpolating u field"
    u_init = interpolate_per_level(u_init, old_degree, new_degree, (Face, Center, Center), H, H_orig)
    @info "interpolating v field"
    v_init = interpolate_per_level(v_init, old_degree, new_degree, (Center, Face, Center), H, H_orig)
    @info "interpolating v field"
    w_init = interpolate_per_level(w_init, old_degree, new_degree, (Center, Center, Face), H, H_orig)
    @info "interpolating η field"
    η_init = interpolate_per_level(η_init, old_degree, new_degree, (Center, Center, Nothing), H, H_orig)
end
