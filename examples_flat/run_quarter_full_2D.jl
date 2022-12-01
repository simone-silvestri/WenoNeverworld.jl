using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb

output_dir    = joinpath(@__DIR__, "../files_four_full_2D")
@show output_prefix = output_dir * "/neverworld_quarter_full_2D"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = false
init_file   = "files_four_isopycnal_new_bathy/neverworld_quarter_isopycnal_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 2minutes
stop_time = 100years

include("../src/new_divergence.jl")

import Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

using Oceananigans.Operators
using Oceananigans.Advection: WENOVectorInvariant, 
		    _left_biased_interpolate_yᵃᶜᵃ,
		    _left_biased_interpolate_xᶜᵃᵃ,
		   _right_biased_interpolate_yᵃᶜᵃ,
		   _right_biased_interpolate_xᶜᵃᵃ,
		    	    upwind_biased_product,
			    VorticityStencil

@inline function vertical_vorticity_U(i, j, k, grid, scheme::WENOVectorInvariant{N, FT, XT, YT, ZT, VI}, u, v) where {N, FT, XT, YT, ZT, VI}
    v̂  =  ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ζᴸ = _center_interpolate_xᶠᵃᵃ(i, j, k, grid,  _left_biased_interpolate_yᵃᶜᵃ, scheme, ζ₃ᶠᶠᶜ, VI, u, v)
    ζᴿ = _center_interpolate_xᶠᵃᵃ(i, j, k, grid, _right_biased_interpolate_yᵃᶜᵃ, scheme, ζ₃ᶠᶠᶜ, VI, u, v)
    return - upwind_biased_product(v̂, ζᴸ, ζᴿ) 
end

@inline function vertical_vorticity_V(i, j, k, grid, scheme::WENOVectorInvariant{N, FT, XT, YT, ZT, VI}, u, v) where {N, FT, XT, YT, ZT, VI}
    û  =  ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ζᴸ = _center_interpolate_yᵃᶠᵃ(i, j, k, grid,  _left_biased_interpolate_xᶜᵃᵃ, scheme, ζ₃ᶠᶠᶜ, VI, u, v)
    ζᴿ = _center_interpolate_yᵃᶠᵃ(i, j, k, grid, _right_biased_interpolate_xᶜᵃᵃ, scheme, ζ₃ᶠᶠᶜ, VI, u, v)
    return + upwind_biased_product(û, ζᴸ, ζᴿ) 
end

tracer_advection = WENO(grid.underlying_grid) 

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, tracer_advection)

increase_simulation_Δt!(simulation, cutoff_time = 75years, new_Δt = 5.0minutes)
increase_simulation_Δt!(simulation, cutoff_time = 76years, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 77years, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = false
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
