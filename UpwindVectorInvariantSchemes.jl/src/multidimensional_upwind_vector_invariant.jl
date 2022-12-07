using Oceananigans.Advection: _advective_momentum_flux_Wu, _advective_momentum_flux_Wv

struct MultiDimensionalVectorInvariant{N, FT, Z, D, V} <: AbstractAdvectionScheme{N, FT}
    ζ_upwind_scheme :: Z
    δ_upwind_scheme :: D
    vertical_scheme :: V

    function MultiDimensionalVectorInvariant{N, FT}(ζ_upwind_scheme::Z, δ_upwind_scheme::D, vertical_scheme::V) where {N, FT, Z, D, V}
        return new{N, FT, Z, D, V}(ζ_upwind_scheme, δ_upwind_scheme, vertical_scheme)
    end
end

validate_momentum_advection(momentum_advection::MultiDimensionalVectorInvariant, grid::AbstractHorizontallyCurvilinearGrid) = momentum_advection

function MultiDimensionalVectorInvariant(; ζ_upwind_scheme::AbstractAdvectionScheme{N, FT} = WENO(VelocityStencil()), 
                                           δ_upwind_scheme = WENO(VorticityStencil()), 
                                           vertical_scheme = Centered(), kwargs...) where {N, FT}

    return MultiDimensionalVectorInvariant{N, FT}(ζ_upwind_scheme, δ_upwind_scheme, vertical_scheme)
end

Adapt.adapt_structure(to, scheme::MultiDimensionalVectorInvariant{N, FT}) where {N, FT} =
MultiDimensionalVectorInvariant{N, FT}(Adapt.adapt(to, scheme.ζ_upwind_scheme), 
                                       Adapt.adapt(to, scheme.δ_upwind_scheme), 
                                       Adapt.adapt(to, scheme.vertical_scheme))

@inline vertical_scheme(scheme::MultiDimensionalVectorInvariant) = string(nameof(typeof(scheme.vertical_scheme)))
@inline smoothness_stencil(::MultiDimensionalVectorInvariant{<:Any, <:Any, <:WENO{N, FT, XT, YT, ZT, VI}}) where {N, FT, XT, YT, ZT, VI} = VI

@inline U_dot_∇u(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, U) = (
    + upwinded_vector_invariant_U(i, j, k, grid, scheme, U.u, U.v)      
    + vertical_advection_U(i, j, k, grid, scheme, U)
    + bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v))              
    
@inline U_dot_∇v(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, U) = (
    + upwinded_vector_invariant_V(i, j, k, grid, scheme, U.u, U.v)      
    + vertical_advection_V(i, j, k, grid, scheme, U)
    + bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v))

@inline bernoulli_head_U(i, j, k, grid, ::MultiDimensionalVectorInvariant, u, v) = ∂xᶠᶜᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)
@inline bernoulli_head_V(i, j, k, grid, ::MultiDimensionalVectorInvariant, u, v) = ∂yᶜᶠᶜ(i, j, k, grid, Khᶜᶜᶜ, u, v)
    
@inline vertical_advection_U(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, U) = 
    1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme.vertical_scheme, U.w, U.u)

@inline vertical_advection_V(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, U) = 
    1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme.vertical_scheme, U.w, U.v)
         
@inline function upwinded_vector_invariant_U(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, u, v)
    
    Sζ = smoothness_stencil(scheme.ζ_upwind_scheme)

    @inbounds v̂ = ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ζᴸ = multi_dimensional_reconstruction_x(i, j, k, grid, scheme.ζ_upwind_scheme,  _left_biased_interpolate_yᵃᶜᵃ, ζ₃ᶠᶠᶜ, Sζ, u, v)
    ζᴿ = multi_dimensional_reconstruction_x(i, j, k, grid, scheme.ζ_upwind_scheme, _right_biased_interpolate_yᵃᶜᵃ, ζ₃ᶠᶠᶜ, Sζ, u, v)

    Sδ = smoothness_stencil(scheme.δ_upwind_scheme)
    
    @inbounds û = u[i, j, k]
    δᴸ = multi_dimensional_reconstruction_y(i, j, k, grid, scheme.δ_upwind_scheme,  _left_biased_interpolate_xᶠᵃᵃ, div_xyᶜᶜᶜ, Sδ, u, v)
    δᴿ = multi_dimensional_reconstruction_y(i, j, k, grid, scheme.δ_upwind_scheme, _right_biased_interpolate_xᶠᵃᵃ, div_xyᶜᶜᶜ, Sδ, u, v)

    return upwind_biased_product(û, δᴸ, δᴿ) - upwind_biased_product(v̂, ζᴸ, ζᴿ)
end

@inline function upwinded_vector_invariant_V(i, j, k, grid, scheme::MultiDimensionalVectorInvariant, u, v) 

    Sζ = smoothness_stencil(scheme.ζ_upwind_scheme)

    @inbounds û  =  ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ζᴸ = multi_dimensional_reconstruction_y(i, j, k, grid, scheme.ζ_upwind_scheme,  _left_biased_interpolate_xᶜᵃᵃ, ζ₃ᶠᶠᶜ, Sζ, u, v)
    ζᴿ = multi_dimensional_reconstruction_y(i, j, k, grid, scheme.ζ_upwind_scheme, _right_biased_interpolate_xᶜᵃᵃ, ζ₃ᶠᶠᶜ, Sζ, u, v)

    Sδ = smoothness_stencil(scheme.δ_upwind_scheme)

    @inbounds v̂ = v[i, j, k]
    δᴸ = multi_dimensional_reconstruction_x(i, j, k, grid, scheme.δ_upwind_scheme,  _left_biased_interpolate_yᵃᶠᵃ, div_xyᶜᶜᶜ, Sδ, u, v)
    δᴿ = multi_dimensional_reconstruction_x(i, j, k, grid, scheme.δ_upwind_scheme, _right_biased_interpolate_yᵃᶠᵃ, div_xyᶜᶜᶜ, Sδ, u, v)

    return upwind_biased_product(û, ζᴸ, ζᴿ) + upwind_biased_product(v̂, δᴸ, δᴿ)
end