struct CenterVerticalScheme end
struct UpwindVerticalScheme end

## Standard vertical advection divided into δ + ∂z(uw)

@inline function vertical_advection_U(i, j, k, grid, scheme::UpwindVectorInvariant, U)

    VI = smoothness_stencil(scheme)

    u, v, W = U

    @inbounds û = u[i, j, k] 
    δᴸ    =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δᴿ    = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δflux = upwind_biased_product(û, δᴸ, δᴿ) 
    wterm =  1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme.upwind_scheme, W, u)

    return wterm + δflux
end

@inline function vertical_advection_V(i, j, k, grid, scheme::UpwindVectorInvariant, U)

    VI = smoothness_stencil(scheme)

    u, v, W = U

    @inbounds v̂ = v[i, j, k] 
    δᴸ    =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δᴿ    = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δflux = upwind_biased_product(v̂, δᴸ, δᴿ) 
    wterm =  1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme.upwind_scheme, W, v)

    return wterm + δflux
end

## Conservative vertical advection

const ConservativeVerticalVectorInvariant = UpwindVectorInvariant{<:Any, <:Any, <:Any, <:CenterVerticalScheme}

@inline vertical_advection_U(i, j, k, grid, scheme::ConservativeVerticalVectorInvariant, U) = vertical_advection_U(i, j, k, grid, scheme.vertical_scheme, U)
@inline vertical_advection_V(i, j, k, grid, scheme::ConservativeVerticalVectorInvariant, U) = vertical_advection_V(i, j, k, grid, scheme.vertical_scheme, U)

@inline vertical_advection_U(i, j, k, grid, ::CenterVerticalScheme, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₂wᶠᶜᶠ, U.u, U.w) / Azᶠᶜᶜ(i, j, k, grid)
@inline vertical_advection_V(i, j, k, grid, ::CenterVerticalScheme, U) =  ℑzᵃᵃᶜ(i, j, k, grid, ζ₁wᶜᶠᶠ, U.v, U.w) / Azᶜᶠᶜ(i, j, k, grid)

@inbounds ζ₂wᶠᶜᶠ(i, j, k, grid, u, w) = ℑxᶠᵃᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶠᶜᶠ(i, j, k, grid, u) 
@inbounds ζ₁wᶜᶠᶠ(i, j, k, grid, v, w) = ℑyᵃᶠᵃ(i, j, k, grid, Az_qᶜᶜᶠ, w) * ∂zᶜᶠᶠ(i, j, k, grid, v) 
    