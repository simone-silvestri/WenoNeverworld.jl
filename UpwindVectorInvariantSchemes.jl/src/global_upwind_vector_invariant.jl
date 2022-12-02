using Oceananigans.Advection: _advective_momentum_flux_Wu, _advective_momentum_flux_Wv

struct GlobalVectorInvariant{N, FT, U, V} <: AbstractAdvectionScheme{N, FT}
    upwind_scheme   :: U
    vertical_scheme :: V

    function GlobalVectorInvariant{N, FT}(upwind_scheme::U, vertical_scheme::V) where {N, FT, U, V}
        return new{N, FT, U, V}(upwind_scheme, vertical_scheme)
    end
end

validate_momentum_advection(momentum_advection::GlobalVectorInvariant, grid::AbstractHorizontallyCurvilinearGrid) = momentum_advection

function GlobalVectorInvariant(; upwind_scheme::AbstractAdvectionScheme{N, FT} = WENO(VelocityStencil()), 
                                 vertical_scheme = Centered(), kwargs...) where {N, FT}

    return GlobalVectorInvariant{N, FT}(upwind_scheme, vertical_scheme)
end

Adapt.adapt_structure(to, scheme::GlobalVectorInvariant{N, FT}) where {N, FT} =
        GlobalVectorInvariant{N, FT}(Adapt.adapt(to, scheme.upwind_scheme), Adapt.adapt(to, vertical_scheme))

@inline vertical_scheme(scheme::GlobalVectorInvariant)       = string(nameof(typeof(scheme.vertical_scheme)))
@inline smoothness_stencil(::GlobalVectorInvariant{<:Any, <:Any, <:WENO{N, FT, XT, YT, ZT, VI}}) where {N, FT, XT, YT, ZT, VI} = VI

@inline U_dot_∇u(i, j, k, grid, scheme::GlobalVectorInvariant, U) = (
    + upwinded_vector_invariant_U(i, j, k, grid, scheme, U.u, U.v)      
    + vertical_advection_U(i, j, k, grid, scheme, U))              
    
@inline U_dot_∇v(i, j, k, grid, scheme::GlobalVectorInvariant, U) = (
    + upwinded_vector_invariant_V(i, j, k, grid, scheme, U.u, U.v)      
    + vertical_advection_V(i, j, k, grid, scheme, U))

@inline vertical_advection_U(i, j, k, grid, scheme::GlobalVectorInvariant, U) = 
    1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme.vertical_scheme, U.w, U.u)

@inline vertical_advection_V(i, j, k, grid, scheme::GlobalVectorInvariant, U) = 
    1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme.vertical_scheme, U.w, U.v)

@inline δ_plus_∂xu(i, j, k, grid, u, v) = div_xyᶜᶜᶜ(i, j, k, grid, u, v) + δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline ζ_plus_∂yu(i, j, k, grid, u, v) =   + ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) + δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, u) / Azᶠᶠᶜ(i, j, k, grid)
@inline δ_plus_∂yv(i, j, k, grid, u, v) = div_xyᶜᶜᶜ(i, j, k, grid, u, v) + δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)
@inline ζ_plus_∂xv(i, j, k, grid, u, v) =   - ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) + δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, v) / Azᶠᶠᶜ(i, j, k, grid)
                    
@inline function upwinded_vector_invariant_U(i, j, k, grid, scheme::GlobalVectorInvariant, u, v)
    
    VI = smoothness_stencil(scheme)

    @inbounds v̂ = ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ζVᴸ =  _left_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme.upwind_scheme, ζ_plus_∂xv, VI, u, v)
    ζVᴿ = _right_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme.upwind_scheme, ζ_plus_∂xv, VI, u, v)

    @inbounds û = u[i, j, k]
    δUᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, δ_plus_∂xu, VI, u, v)
    δUᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, δ_plus_∂xu, VI, u, v)

    return upwind_biased_product(û, δUᴸ, δUᴿ) + upwind_biased_product(v̂, ζVᴸ, ζVᴿ)
end

@inline function upwinded_vector_invariant_V(i, j, k, grid, scheme::GlobalVectorInvariant, u, v) 

    VI = smoothness_stencil(scheme)

    @inbounds û  =  ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ζUᴸ =  _left_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ζ_plus_∂yu, VI, u, v)
    ζUᴿ = _right_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ζ_plus_∂yu, VI, u, v)

    @inbounds v̂ = v[i, j, k]
    δVᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, δ_plus_∂yv, VI, u, v)
    δVᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, δ_plus_∂yv, VI, u, v)

    return upwind_biased_product(û, ζUᴸ, ζUᴿ) + upwind_biased_product(v̂, δVᴸ, δVᴿ)
end