using Oceananigans
using Oceananigans.Grids: halo_size, xnode, ynode, znode
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.Advection: VelocityStencil
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb
using Statistics: mean, quantile
using JLD2

output_dir    = joinpath(@__DIR__, "../files_spin_down")

arch   = GPU()
new_degree = 1/4

grid     = NeverworldGrid(arch,  1/4)
grid_aux = NeverworldGrid(CPU(), 1/4)
Hx, Hy, Hz = halo_size(grid_aux)

# Extend the vertical advection scheme
interp_init = false
init_file   = "files_four_new_bathy/neverworld_quarter_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 10minutes
stop_time = 100days

include("../src/isopycnally_rotated_upwinding.jl")

tadv2 = WENO()
tadv4 = IsopycnallyRotatedUpwindScheme(tadv2, Centered(order = 6))

tracer_advections = [tadv2, tadv4]

momentum_advections = [VectorInvariant(), WENO(vector_invariant = VelocityStencil())]

@inline initialize_tracer(x, y, z) = y > - 60 && y < - 10 && x > 20 && x < 50 && z > - 3000 && z < - 200

using WenoNeverworld: geometric_νhb, default_biharmonic_viscosity

vertical_diffusivity = VerticalScalarDiffusivity(ν = 1e-4, κ = (; b = 1e-5, c = 1e-5))

for (idx_mom, momentum_advection) in enumerate(momentum_advections), (idx_trac, tracer_advection) in enumerate(tracer_advections)

    if idx_mom == 1
        biharmonic_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=geometric_νhb, discrete_form=true, parameters = 5days)
    else
        biharmonic_viscosity = default_biharmonic_viscosity
    end

    # Construct the neverworld simulation
    simulation = weno_neverworld_simulation(; grid, Δt, stop_time, init_file, tracer_advection, momentum_advection,
                                            tracers = (:b, :c), vertical_diffusivity, biharmonic_viscosity)

    b_init = jldopen(init_file)["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    η_init = jldopen(init_file)["η/data"][Hx+1:end-Hx, Hy+1:end-Hy, :]

    Gb⁻_init = jldopen(init_file)["timestepper/G⁻/b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    Gu⁻_init = jldopen(init_file)["timestepper/G⁻/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    Gv⁻_init = jldopen(init_file)["timestepper/G⁻/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]

    Gbⁿ_init = jldopen(init_file)["timestepper/G⁻/b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    Guⁿ_init = jldopen(init_file)["timestepper/G⁻/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    Gvⁿ_init = jldopen(init_file)["timestepper/G⁻/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]

    b_mean = quantile(b_init[:], 0.7)
    c_init = zeros(size(b_init)...)

    for i in 1:size(grid_aux, 1), j in 1:size(grid_aux, 2), k in 1:size(grid_aux, 3)

        x = xnode(Center(), Center(), Center(), i, j, k, grid_aux)
        y = ynode(Center(), Center(), Center(), i, j, k, grid_aux)
        z = znode(Center(), Center(), Center(), i, j, k, grid_aux)

        if initialize_tracer(x, y, z)
            c_init[i, j, k] = exp( - (b_init[i, j, k] - b_mean)^2 / 0.00002)
        end
    end

    model = simulation.model
    
    set!(model.tracers.b, b_init)
    set!(model.tracers.c, c_init)

    set!(model.velocities.u, u_init)
    set!(model.velocities.v, v_init)
    set!(model.velocities.w, w_init)

    set!(model.free_surface.η, η_init)

    set!(model.timestepper.G⁻.b, Gb⁻_init)
    set!(model.timestepper.G⁻.u, Gu⁻_init)
    set!(model.timestepper.G⁻.v, Gv⁻_init)

    set!(model.timestepper.Gⁿ.b, Gbⁿ_init)
    set!(model.timestepper.Gⁿ.u, Guⁿ_init)
    set!(model.timestepper.Gⁿ.v, Gvⁿ_init)

    # Let's goo!
    @info "Running with Δt = $(prettytime(simulation.Δt))"
    @show output_prefix = output_dir * "/spin_down_trac$(idx_trac)_mom$(idx_mom)"

    checkpoint_outputs!(simulation, output_prefix; checkpoint_time = 10days)

    # initializing the time for wall_time calculation
    run_simulation!(simulation; interp_init, init_file = nothing)
end

