struct UpwindKineticScheme end

const UpwindKineticVectorInvariant = UpwindVectorInvariant{<:Any, <:Any, <:Any, <:Any, UpwindKineticScheme}

@inline ∂x_u_c(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u)
@inline ∂x_v_f(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v)
@inline ∂y_u_f(i, j, k, grid, u, v) = ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline ∂y_v_c(i, j, k, grid, u, v) = ∂yᶜᶜᶜ(i, j, k, grid, v)

@inline function bernoulli_head_U(i, j, k, grid, scheme::UpwindKineticVectorInvariant, u, v) 

    VI = smoothness_stencil(scheme)

    @inbounds û = u[i, j, k]
    ∂uᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, ∂x_u_c, VI, u, v)
    ∂uᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, ∂x_u_c, VI, u, v)

    @inbounds v̂ = ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ∂vᴸ =  _left_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme, ∂x_v_f, VI, u, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme, ∂x_v_f, VI, u, v)

    return upwind_biased_product(û, ∂uᴸ, ∂uᴿ) + upwind_biased_product(v̂, ∂vᴸ, ∂vᴿ)
end

@inline function bernoulli_head_V(i, j, k, grid, scheme::UpwindKineticVectorInvariant, u, v) 

    VI = smoothness_stencil(scheme)

    @inbounds û = ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ∂uᴸ =  _left_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme, ∂y_u_f, VI, u, v)
    ∂uᴿ = _right_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme, ∂y_u_f, VI, u, v)

    @inbounds v̂ = v[i, j, k]
    ∂vᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, ∂y_v_c, VI, u, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, ∂y_v_c, VI, u, v)

    return upwind_biased_product(û, ∂uᴸ, ∂uᴿ) + upwind_biased_product(v̂, ∂vᴸ, ∂vᴿ)
end