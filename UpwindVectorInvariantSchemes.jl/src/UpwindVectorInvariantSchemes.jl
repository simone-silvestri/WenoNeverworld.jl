module UpwindVectorInvariantSchemes

export UpwindVectorInvariant, GlobalVectorInvariant, UpwindKineticScheme, UpwindVerticalScheme

using Oceananigans
using WenoNeverworld
using Oceananigans.Advection: AbstractAdvectionScheme
using Oceananigans.Advection: VorticityStencil

import Oceananigans.Advection: vertical_advection_U, vertical_advection_V
import Oceananigans.Advection: bernoulli_head_U, bernoulli_head_V
import Oceananigans.Advection: U_dot_∇u, U_dot_∇v
import Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

@inline smoothness_stencil(scheme) = VorticityStencil

import Oceananigans.Advection: WENO
using Oceananigans.Advection
using Oceananigans.Advection: compute_reconstruction_coefficients, SmoothnessStencil

using Oceananigans.Operators

function WENO(vector_invariant::SmoothnessStencil, FT::DataType=Float64; 
              order = 5,
              grid = nothing, 
              zweno = true,
              bounds = nothing)

    if !(grid isa Nothing) 
        FT = eltype(grid)
    end

    mod(order, 2) == 0 && throw(ArgumentError("WENO reconstruction scheme is defined only for odd orders"))

    if order < 3
        # WENO(order = 1) is equivalent to UpwindBiased(order = 1)
        return UpwindBiased(order = 1)
    else
        VI = typeof(vector_invariant)
        N  = Int((order + 1) ÷ 2)

        weno_coefficients = compute_reconstruction_coefficients(grid, FT, :WENO; order = N)
        buffer_scheme     = WENO(vector_invariant, FT; grid, order = order - 2, zweno, bounds)
        advecting_velocity_scheme = Centered(FT; grid, order = order - 1)
    end

    return WENO{N, FT, VI, zweno}(weno_coefficients..., bounds, buffer_scheme, advecting_velocity_scheme)
end

include("default_vector_invariant.jl")
include("upwinded_kinetic_energy_scheme.jl")
include("upwinded_vertical_advection.jl")
include("global_upwind_vector_invariant.jl")

end 