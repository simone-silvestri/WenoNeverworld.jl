struct UpwindVerticalScheme end

const UpwindVerticalVectorInvariant = UpwindVectorInvariant{<:Any, <:Any, <:Any, UpwindVerticalScheme}

@inline function vertical_advection_U(i, j, k, grid, scheme::UpwindVerticalVectorInvariant, U)

    VI = smoothness_stencil(scheme)

    u, v, W = U

    @inbounds û = u[i, j, k] 
    δᴸ    =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δᴿ    = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δflux = upwind_biased_product(û, δᴸ, δᴿ) 
    wterm =  1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme.upwind_scheme, W, u)

    return wterm + δflux
end

@inline function vertical_advection_V(i, j, k, grid, scheme::UpwindVerticalVectorInvariant, U)

    VI = smoothness_stencil(scheme)

    u, v, W = U

    @inbounds v̂ = v[i, j, k] 
    δᴸ    =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δᴿ    = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, div_xyᶜᶜᶜ, VI, u, v)
    δflux = upwind_biased_product(v̂, δᴸ, δᴿ) 
    wterm =  1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme.upwind_scheme, W, v)

    return wterm + δflux
end
