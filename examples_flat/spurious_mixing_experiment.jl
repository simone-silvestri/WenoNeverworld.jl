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

grid = NeverworldGrid(arch, 1/4)
Hx, Hy, Hz = halo_size(grid)

# Extend the vertical advection scheme
interp_init = false
init_file   = "files_four_new_bathy/neverworld_quarter_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 10minutes
stop_time = 100days

include("../src/isopycnally_rotated_upwinding.jl")

tadv1 = WENO(grid.underlying_grid)
tadv2 = WENO()
tadv3 = IsopycnallyRotatedUpwindScheme(tadv1, Centered(grid.underlying_grid, order = 6))
tadv4 = IsopycnallyRotatedUpwindScheme(tadv2, Centered(order = 6))

tracer_advections = [tadv1, tadv2, tadv3, tadv4]

momentum_advections = [VectorInvariant(), WENO(vector_invariant = VelocityStencil())]

@inline initialize_tracer(x, y, z) = y > - 60 && y < - 10 && x > 30 && x < 50 && z > - 3000 && z < - 200

using WenoNeverworld: geometric_νhb, default_biharmonic_viscosity

vertical_diffusivity = VerticalScalarDiffusivity(ν = 1e-4, κ = (; b = 1e-5, c = 0.0))

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
    b_mean = quantile(b_init[:], 0.7)
    c_init = zeros(size(b_init)...)

    for i in 1:size(grid, 1), j in 1:size(grid, 2), k in 1:size(grid, 3)

        x = xnode(Center(), Center(), Center(), i, j, k, grid)
        y = ynode(Center(), Center(), Center(), i, j, k, grid)
        z = znode(Center(), Center(), Center(), i, j, k, grid)

        if initialize_tracer(x, y, z)
            c_init[i, j, k] = exp( - (b_init[i, j, k] - b_mean)^2 / 0.00002)
        end
    end
                                            
    set!(simulation.model.tracers.c, c_init)
    # Let's goo!
    @info "Running with Δt = $(prettytime(simulation.Δt))"
    @show output_prefix = output_dir * "/spin_down_trac$(idx_trac)_mom$(idx_mom)"

    checkpoint_outputs!(simulation, output_prefix; checkpoint_time = 10days)

    # initializing the time for wall_time calculation
    run_simulation!(simulation; interp_init, init_file)
end

