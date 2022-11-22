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

using Oceananigans.Operators
using Oceananigans.Operators: Vᶠᶜᶜ, Vᶜᶠᶜ, δzᵃᵃᶜ, Azᶠᶜᶠ, Azᶜᶠᶠ, ℑxzᶠᵃᶜ, ℑyzᵃᶠᶜ

import Oceananigans.Advection: vertical_advection_U, vertical_advection_V, advective_momentum_flux_Wu, advective_momentum_flux_Wv
    
@inline vertical_advection_U(i, j, k, grid, ::WENOVectorInvariant, u, w) =  ℑxzᶠᵃᶜ(i, j, k, grid, w) * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)
@inline vertical_advection_V(i, j, k, grid, ::WENOVectorInvariant, v, w) =  ℑyzᵃᶠᶜ(i, j, k, grid, w) * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v)

interp_init = true
init_file   = "files_lowres_new_bathy/neverworld_lowres_checkpoint_iteration2067840.jld2" 

# Simulation parameters
Δt        = 1minutes
stop_time = 7000days

biharmonic_viscosity = HorizontalDivergenceScalarDiffusivity(ν = 100.0)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, biharmonic_viscosity, vertical_diffusivity)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 2.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 100days, new_Δt = 5.0minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 300days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

standard_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
