using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity

output_dir    = joinpath(@__DIR__, "../files_four_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
using Oceananigans.Advection: WENOVectorInvariant, _advective_momentum_flux_Wu, _advective_momentum_flux_Wv
using Oceananigans.Advection:  _left_biased_interpolate_xᶠᵃᵃ, 
                              _right_biased_interpolate_xᶠᵃᵃ,
                               _left_biased_interpolate_yᵃᶠᵃ,
                              _right_biased_interpolate_yᵃᶠᵃ,
                               _left_biased_interpolate_zᵃᵃᶠ,
                              _right_biased_interpolate_zᵃᵃᶠ,
                                        upwind_biased_product

using Oceananigans.Operators: Vᶠᶜᶜ, Vᶜᶠᶜ, δzᵃᵃᶜ, Azᶠᶜᶠ, Azᶜᶠᶠ

import Oceananigans.Advection: vertical_advection_U, vertical_advection_V, advective_momentum_flux_Wu, advective_momentum_flux_Wv

@inline function advective_momentum_flux_Wu(i, j, k, grid, scheme::WENOVectorInvariant, W, u)

    wᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, W)
    wᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, W)
    uᴸ =  _left_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, u)
    uᴿ = _right_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, u)

    w̃ = 0.5 * (wᴸ + wᴿ)

    return Azᶠᶜᶠ(i, j, k, grid) * upwind_biased_product(w̃, uᴸ, uᴿ)
end

@inline function advective_momentum_flux_Wv(i, j, k, grid, scheme::WENOVectorInvariant, W, v)

    wᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, W)
    wᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, W)
    vᴸ =  _left_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, v)
    vᴿ = _right_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme, v)

    w̃ = 0.5 * (wᴸ + wᴿ)

    return Azᶜᶠᶠ(i, j, k, grid) * upwind_biased_product(w̃, vᴸ, vᴿ)
end

@inline vertical_advection_U(i, j, k, grid, scheme::WENOVectorInvariant, u, w) = 1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme, w, u)
@inline vertical_advection_V(i, j, k, grid, scheme::WENOVectorInvariant, v, w) = 1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme, w, v)

interp_init = true
init_file   = "files_lowres_new_bathy/neverworld_lowres_checkpoint_iteration2067840.jld2" 

# Simulation parameters
Δt        = 2minutes
stop_time = 7000days

biharmonic_viscosity = HorizontalDivergenceScalarDiffusivity(ν = 100.0)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, biharmonic_viscosity, vertical_diffusivity)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 100days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
