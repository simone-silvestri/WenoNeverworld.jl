using Oceananigans
using Oceananigans.Operators
using Oceananigans.Advection: UpwindScheme

using Oceananigans.Advection: _left_biased_interpolate_xᶠᵃᵃ,
                             _right_biased_interpolate_xᶠᵃᵃ,
                              _left_biased_interpolate_yᵃᶠᵃ,
                             _right_biased_interpolate_yᵃᶠᵃ,
                                      upwind_biased_product

import Oceananigans.Models.HydrostaticFreeSurfaceModels: flux_div_xyᶜᶜᶠ
import Oceananigans.Operators: div_xyᶜᶜᶜ, div_xyᶜᶜᶠ
import Oceananigans.Advection: advective_tracer_flux_x, advective_tracer_flux_y, advective_tracer_flux_z

const ε = 1e-6
const two_32 = Int32(2)

"""
    Third order center to center reconstruction coefficients (left biased, centered, right biased) obtained with

`stencil_coefficients(50, 0, collect(0.5:99.5), collect(1:100); order = 3)`
`stencil_coefficients(50, 1, collect(0.5:99.5), collect(1:100); order = 3)`
`stencil_coefficients(50, 2, collect(0.5:99.5), collect(1:100); order = 3)`
"""
const C₀ = (-1.0, 2.0,  23.0) ./ 24
const C₁ = (-1.0, 26.0, -1.0) ./ 24
const C₂ = (23.0, 2.0,  -1.0) ./ 24

"""
    Fifth order center to center reconstruction coefficients (centered) obtained with

`stencil_coefficients(50, 2, collect(0.5:99.5), collect(1:100); order = 5)`
"""
fifth_order_coeffs = (9.0, -116.0, 2134.0, -116.0, 9.0) ./ 1920

""" 
Optimal coefficients obtained assuming all stencils are perfectly smooth (`β₀ = β₁ = β₂`)
and the reconstruction is fifth order
"""
const σ⁺ = 214/80
const σ⁻ =  67/40

const OC₀⁺ =   9.0 / 80 / σ⁺
const OC₁⁺ =  49.0 / 20 / σ⁺
const OC₂⁺ =   9.0 / 80 / σ⁺

const OC₀⁻ =   9.0 / 40 / σ⁻
const OC₁⁻ =  49.0 / 40 / σ⁻
const OC₂⁻ =   9.0 / 40 / σ⁻

"""
To calculate the smoothness coefficients

    `βᵣ = ∑ₗ₌₁² ∫ₓ₋ˣ⁺ Δx²ˡ⁻¹ ∂ˡpᵣ(x) dx`

where

`pᵣ(x) = ∑ₘ₌₀³∑ⱼ₌₀ᵐ⁻¹ v̂ᵢ₋ᵣ₊ⱼ Δx`
"""
@inline   left_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * ( ψ[1] - 4ψ[2] + 3ψ[3])^two_32
@inline center_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * ( ψ[1]         -  ψ[3])^two_32
@inline  right_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * (3ψ[1] - 4ψ[2] +  ψ[3])^two_32

@inline function centered_reconstruction_weights(FT, ψ₀, ψ₁, ψ₂)

    β₀ =   left_biased_β(FT, ψ₀)
    β₁ = center_biased_β(FT, ψ₁)
    β₂ =  right_biased_β(FT, ψ₂)

    α₀⁺ = FT(OC₀⁺) / (β₀ + FT(ε))^two_32
    α₁⁺ = FT(OC₁⁺) / (β₁ + FT(ε))^two_32
    α₂⁺ = FT(OC₂⁺) / (β₂ + FT(ε))^two_32

    α₀⁻ = FT(OC₀⁻) / (β₀ + FT(ε))^two_32
    α₁⁻ = FT(OC₁⁻) / (β₁ + FT(ε))^two_32
    α₂⁻ = FT(OC₂⁻) / (β₂ + FT(ε))^two_32

    Σα⁺ = α₀⁺ + α₁⁺ + α₂⁺
    w₀⁺ = α₀⁺ / Σα⁺
    w₁⁺ = α₁⁺ / Σα⁺
    w₂⁺ = α₂⁺ / Σα⁺

    Σα⁻ = α₀⁻ + α₁⁻ + α₂⁻
    w₀⁻ = α₀⁻ / Σα⁻
    w₁⁻ = α₁⁻ / Σα⁻
    w₂⁻ = α₂⁻ / Σα⁻

    return w₀⁺, w₁⁺, w₂⁺, w₀⁻, w₁⁻, w₂⁻
end

@inline function weno_reconstruction(FT, ψ₀, ψ₁, ψ₂)

    w₀⁺, w₁⁺, w₂⁺, w₀⁻, w₁⁻, w₂⁻ = centered_reconstruction_weights(FT, ψ₀, ψ₁, ψ₂)

    ψ̂₀ = sum(ψ₀ .* C₀)
    ψ̂₁ = sum(ψ₁ .* C₁)
    ψ̂₂ = sum(ψ₂ .* C₂)

    return σ⁺ * (ψ̂₀ * w₀⁺ + ψ̂₁ * w₁⁺ + ψ̂₂ * w₂⁺) - 
           σ⁻ * (ψ̂₀ * w₀⁻ + ψ̂₁ * w₁⁻ + ψ̂₂ * w₂⁻) 
end

@inline function center_interpolate_xᶠᵃᵃ(i, j, k, grid, u)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (u[i-2, j, k], u[i-1, j, k], u[i,   j, k])
    @inbounds ψ₁ = (u[i-1, j, k], u[i,   j, k], u[i+1, j, k])
    @inbounds ψ₂ = (u[i,   j, k], u[i+1, j, k], u[i+2, j, k])
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂)
end

@inline function center_interpolate_xᶠᵃᵃ(i, j, k, grid, f::Function, args...)

    FT = eltype(grid)

    @inbounds ψ₀ = (f(i-2, j, k, grid, args...), f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...))
    @inbounds ψ₁ = (f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...))
    @inbounds ψ₂ = (f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...), f(i+2, j, k, grid, args...))
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂)
end

@inline function center_interpolate_yᵃᶠᵃ(i, j, k, grid, v)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (v[i, j-2, k], v[i, j-1, k], v[i, j,   k])
    @inbounds ψ₁ = (v[i, j-1, k], v[i, j,   k], v[i, j+1, k])
    @inbounds ψ₂ = (v[i, j,   k], v[i, j+1, k], v[i, j+2, k])
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂)
end

@inline function center_interpolate_yᵃᶠᵃ(i, j, k, grid, f::Function, args...)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (f(i, j-2, k, grid, args...), f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...))
    @inbounds ψ₁ = (f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...))
    @inbounds ψ₂ = (f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...), f(i, j+2, k, grid, args...))
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂)
end

using Oceananigans.Grids: inactive_node
using Oceananigans.ImmersedBoundaries: calc_inactive_stencil

const f = Face()
const c = Center()

@eval begin
    @inline near_boundary_center_xᶠᵃᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(3, :left, :x, :ᶜ, xside = :ᶜ)...))
    @inline near_boundary_center_yᵃᶠᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(3, :left, :y, :ᶜ, yside = :ᶜ)...))
end

@inline _center_interpolate_yᵃᶠᵃ(i, j, k, grid, v) = 
        ifelse(near_boundary_center_yᵃᶠᵃ(i, j, k, grid),  v[i, j, k], center_interpolate_yᵃᶠᵃ(i, j, k, grid, v))

@inline _center_interpolate_xᶠᵃᵃ(i, j, k, grid, u) = 
        ifelse(near_boundary_center_xᶠᵃᵃ(i, j, k, grid),  u[i, j, k], center_interpolate_xᶠᵃᵃ(i, j, k, grid, u))

@inline _center_interpolate_yᵃᶠᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_yᵃᶠᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_yᵃᶠᵃ(i, j, k, grid, f, args...))

@inline _center_interpolate_xᶠᵃᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_xᶠᵃᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_xᶠᵃᵃ(i, j, k, grid, f, args...))

@inline flux_div_xyᶜᶜᶠ(i, j, k, grid, Qu, Qv) = δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᵃᶠᵃ, Qu) + δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᶠᵃᵃ, Qv)

@inline div_xyᶜᶜᶜ(i, j, k, grid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᵃᶠᵃ, Δy_qᶠᶜᶜ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᶠᵃᵃ, Δx_qᶜᶠᶜ, v))

@inline div_xyᶜᶜᶜ(i, j, k, grid::ImmersedBoundaryGrid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᵃᶠᵃ, Δy_qᶠᶜᶜ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᶠᵃᵃ, Δx_qᶜᶠᶜ, v))

@inline div_xyᶜᶜᶠ(i, j, k, grid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᵃᶠᵃ, Δy_qᶠᶜᶠ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᶠᵃᵃ, Δx_qᶜᶠᶠ, v))

@inline div_xyᶜᶜᶠ(i, j, k, grid::ImmersedBoundaryGrid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᵃᶠᵃ, Δy_qᶠᶜᶠ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᶠᵃᵃ, Δx_qᶜᶠᶠ, v))
                            
@inline function advective_tracer_flux_x(i, j, k, grid, scheme::UpwindScheme, U, c) 

    ũ  =       _center_interpolate_yᵃᶠᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, U)
    cᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, c)
    cᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, c)

    return Axᶠᶜᶜ(i, j, k, grid) * upwind_biased_product(ũ, cᴸ, cᴿ)
end

@inline function advective_tracer_flux_y(i, j, k, grid, scheme::UpwindScheme, V, c)

    ṽ  =       _center_interpolate_xᶠᵃᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, V)
    cᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, c)
    cᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, c)

    return Ayᶜᶠᶜ(i, j, k, grid) * upwind_biased_product(ṽ, cᴸ, cᴿ)
end
