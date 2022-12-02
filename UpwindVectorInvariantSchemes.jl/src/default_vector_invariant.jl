struct CenterVerticalScheme end
struct CenterKineticScheme end

struct UpwindVectorInvariant{N, FT, U, V, K} <: AbstractAdvectionScheme{N, FT}
    upwind_scheme         :: U
    vertical_scheme       :: V

    function UpwindVectorInvariant{N, FT}(upwind_scheme::U, vertical_scheme::V) where {N, FT, U, V}
        return new{N, FT, U, V, K}(upwind_scheme, vertical_scheme)
    end
end

function UpwindVectorInvariant(; upwind_scheme::AbstractAdvectionScheme{N, FT} = VectorInvariant(), 
                                 vertical_scheme = CenterVerticalScheme()) where {N, FT}

    return UpwindVectorInvariant{N, FT}(upwind_scheme, vertical_scheme)
end

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

####
#### Vertical advection terms
####

@inline vertical_advection_U(i, j, k, grid, scheme::ConservativeVerticalVectorInvariant, U) = vertical_advection_U(i, j, k, grid, scheme.vertical_scheme, U)
@inline vertical_advection_V(i, j, k, grid, scheme::ConservativeVerticalVectorInvariant, U) = vertical_advection_V(i, j, k, grid, scheme.vertical_scheme, U)

@inline vertical_advection_U(i, j, k, grid, ::CenterVerticalScheme, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₂wᶠᶜᶠ, U.u, U.w) / Azᶠᶜᶜ(i, j, k, grid)
@inline vertical_advection_V(i, j, k, grid, ::CenterVerticalScheme, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₁wᶜᶠᶠ, U.v, U.w) / Azᶜᶠᶜ(i, j, k, grid)

@inbounds ζ₂wᶠᶜᶠ(i, j, k, grid, u, w) = ℑxᶠᵃᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶠᶜᶠ(i, j, k, grid, u) 
@inbounds ζ₁wᶜᶠᶠ(i, j, k, grid, v, w) = ℑyᵃᶠᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶜᶠᶠ(i, j, k, grid, v) 
    