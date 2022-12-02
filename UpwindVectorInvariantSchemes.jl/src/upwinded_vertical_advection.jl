struct UpwindVerticalScheme end
struct ConservativeUpwindVerticalScheme end

const UpwindVerticalVectorInvariant             = UpwindVectorInvariant{<:Any, <:Any, <:Any, UpwindVerticalScheme}
const ConservativeUpwindVerticalVectorInvariant = UpwindVectorInvariant{<:Any, <:Any, <:Any, ConservativeUpwindVerticalScheme}

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

@inline function bernoulli_head_U(i, j, k, grid, scheme::ConservativeUpwindKineticVectorInvariant, u, v) 

    VI = smoothness_stencil(scheme)

    @inbounds û = u[i, j, k]
    ∂uᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_u_c, VI, u, v)
    ∂uᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_u_c, VI, u, v)
    ∂uˢ =    _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_u_c, VI, u, v)

    @inbounds v̂ = ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ∂vᴸ =  _left_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_v_f, VI, u, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_v_f, VI, u, v)
    ∂vˢ =    _symmetric_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme.upwind_scheme, ∂x_v_f, VI, u, v)

    conservative_term = bernoulli_head_U(i, j, k, grid, CenterKineticScheme(), u, v)

    return conservative_term + 
          (upwind_biased_product(û, ∂uᴸ, ∂uᴿ) - û * ∂uˢ)  + 
          (upwind_biased_product(v̂, ∂vᴸ, ∂vᴿ) - v̂ * ∂vˢ)
end

@inline function bernoulli_head_V(i, j, k, grid, scheme::ConservativeUpwindKineticVectorInvariant, u, v) 

    VI = smoothness_stencil(scheme)

    @inbounds û = ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ∂uᴸ =  _left_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_u_f, VI, u, v)
    ∂uᴿ = _right_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_u_f, VI, u, v)
    ∂uˢ =    _symmetric_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_u_f, VI, u, v)

    @inbounds v̂ = v[i, j, k]
    ∂vᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_v_c, VI, u, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_v_c, VI, u, v)
    ∂vˢ =    _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme, ∂y_v_c, VI, u, v)

    return conservative_term + 
           (upwind_biased_product(û, ∂uᴸ, ∂uᴿ) - û * ∂uˢ)  + 
           (upwind_biased_product(v̂, ∂vᴸ, ∂vᴿ) - v̂ * ∂vˢ)
end