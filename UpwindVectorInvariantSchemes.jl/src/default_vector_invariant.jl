struct DefaultVectorInvariant{N, FT, U} <: AbstractAdvectionScheme{N, FT}
    upwind_scheme         :: U

    function DefaultVectorInvariant{N, FT}(upwind_scheme::U) where {N, FT, U, V}
        return new{N, FT, U, V}(upwind_scheme)
    end
end

using Oceananigans.Grids: AbstractHorizontallyCurvilinearGrid

validate_momentum_advection(momentum_advection::DefaultVectorInvariant, grid::AbstractHorizontallyCurvilinearGrid) = momentum_advection

function DefaultVectorInvariant(; upwind_scheme::AbstractAdvectionScheme{N, FT} = VectorInvariant()) where {N, FT}

    return DefaultVectorInvariant{N, FT}(upwind_scheme)
end

Adapt.adapt_structure(to, scheme::DefaultVectorInvariant{N, FT}) where {N, FT} =
        DefaultVectorInvariant{N, FT}(Adapt.adapt(to, scheme.upwind_scheme), Adapt.adapt(to, vertical_scheme))

@inline smoothness_stencil(::DefaultVectorInvariant{<:Any, <:Any, <:WENO{N, FT, XT, YT, ZT, VI}}) where {N, FT, XT, YT, ZT, VI} = VI

@inline U_dot_∇u(i, j, k, grid, scheme::DefaultVectorInvariant, U) = (
    + vertical_vorticity_U(i, j, k, grid, scheme.upwind_scheme, U.u, U.v)      # Vertical relative vorticity term
    + vertical_advection_U(i, j, k, grid, scheme, U)                               # Horizontal vorticity / vertical advection term
    + bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v))                           # Bernoulli head term
    
@inline U_dot_∇v(i, j, k, grid, scheme::DefaultVectorInvariant, U) = (
    + vertical_vorticity_V(i, j, k, grid, scheme.upwind_scheme, U.u, U.v)      # Vertical relative vorticity term
    + vertical_advection_V(i, j, k, grid, scheme, U)                               # Horizontal vorticity / vertical advection term
    + bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v))                           # Bernoulli head term

@inline bernoulli_head_U(i, j, k, grid, ::DefaultVectorInvariant, u, v) = ∂xᶠᶜᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)
@inline bernoulli_head_V(i, j, k, grid, ::DefaultVectorInvariant, u, v) = ∂yᶜᶠᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)

@inline ϕ²(i, j, k, grid, ϕ)       = @inbounds ϕ[i, j, k]^2
@inline Khᶜᶜᶜ(i, j, k, grid, u, v) = (ℑxᶜᵃᵃ(i, j, k, grid, ϕ², u) + ℑyᵃᶜᵃ(i, j, k, grid, ϕ², v)) / 2

@inline vertical_advection_U(i, j, k, grid, ::DefaultVectorInvariant, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₂wᶠᶜᶠ, U.u, U.w) / Azᶠᶜᶜ(i, j, k, grid)
@inline vertical_advection_V(i, j, k, grid, ::DefaultVectorInvariant, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₁wᶜᶠᶠ, U.v, U.w) / Azᶜᶠᶜ(i, j, k, grid)

@inbounds ζ₂wᶠᶜᶠ(i, j, k, grid, u, w) = ℑxᶠᵃᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶠᶜᶠ(i, j, k, grid, u) 
@inbounds ζ₁wᶜᶠᶠ(i, j, k, grid, v, w) = ℑyᵃᶠᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶜᶠᶠ(i, j, k, grid, v) 
    