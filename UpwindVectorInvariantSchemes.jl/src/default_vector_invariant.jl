struct UpwindVectorInvariant{N, FT, U, V} <: AbstractAdvectionScheme{N, FT}
    upwind_scheme         :: U
    vertical_scheme       :: V

    function UpwindVectorInvariant{N, FT}(upwind_scheme::U, vertical_scheme::V) where {N, FT, U, V}
        return new{N, FT, U, V}(upwind_scheme, vertical_scheme)
    end
end

using Oceananigans.Grids: AbstractHorizontallyCurvilinearGrid

validate_momentum_advection(momentum_advection::UpwindVectorInvariant, grid::AbstractHorizontallyCurvilinearGrid) = momentum_advection

function UpwindVectorInvariant(; upwind_scheme::AbstractAdvectionScheme{N, FT} = VectorInvariant(), 
                                 vertical_scheme = CenterVerticalScheme()) where {N, FT}

    return UpwindVectorInvariant{N, FT}(upwind_scheme, vertical_scheme)
end

Adapt.adapt_structure(to, scheme::UpwindVectorInvariant{N, FT}) where {N, FT} =
        UpwindVectorInvariant{N, FT}(Adapt.adapt(to, scheme.upwind_scheme), Adapt.adapt(to, vertical_scheme))

@inline vertical_scheme(scheme::UpwindVectorInvariant) = string(nameof(typeof(scheme.vertical_scheme)))
@inline smoothness_stencil(::UpwindVectorInvariant{<:Any, <:Any, <:WENO{N, FT, XT, YT, ZT, VI}}) where {N, FT, XT, YT, ZT, VI} = VI

@inline U_dot_∇u(i, j, k, grid, scheme::UpwindVectorInvariant, U) = (
    + vertical_vorticity_U(i, j, k, grid, scheme.upwind_scheme, U.u, U.v)      # Vertical relative vorticity term
    + vertical_advection_U(i, j, k, grid, scheme, U)                               # Horizontal vorticity / vertical advection term
    + bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v))                           # Bernoulli head term
    
@inline U_dot_∇v(i, j, k, grid, scheme::UpwindVectorInvariant, U) = (
    + vertical_vorticity_V(i, j, k, grid, scheme.upwind_scheme, U.u, U.v)      # Vertical relative vorticity term
    + vertical_advection_V(i, j, k, grid, scheme, U)                               # Horizontal vorticity / vertical advection term
    + bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v))                           # Bernoulli head term

@inline bernoulli_head_U(i, j, k, grid, ::UpwindVectorInvariant, u, v) = ∂xᶠᶜᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)
@inline bernoulli_head_V(i, j, k, grid, ::UpwindVectorInvariant, u, v) = ∂yᶜᶠᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)

@inline ϕ²(i, j, k, grid, ϕ)       = @inbounds ϕ[i, j, k]^2
@inline Khᶜᶜᶜ(i, j, k, grid, u, v) = (ℑxᶜᵃᵃ(i, j, k, grid, ϕ², u) + ℑyᵃᶜᵃ(i, j, k, grid, ϕ², v)) / 2
