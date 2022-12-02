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

@inline function centered_reconstruction_weights(FT, s₀, s₁, s₂)

    β₀ =   left_biased_β(FT, s₀)
    β₁ = center_biased_β(FT, s₁)
    β₂ =  right_biased_β(FT, s₂)

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

@inline function weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)

    w₀⁺, w₁⁺, w₂⁺, w₀⁻, w₁⁻, w₂⁻ = centered_reconstruction_weights(FT, s₀, s₁, s₂)

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
    
    @inbounds s₀ = (∂xᶜᶜᶜ(i-2, j, k, grid, ∂xᶜᶜᶜ, u), ∂xᶜᶜᶜ(i-1, j, k, grid, u), ∂xᶜᶜᶜ(i,   j, k, grid, u))
    @inbounds s₁ = (∂xᶜᶜᶜ(i-1, j, k, grid, ∂xᶜᶜᶜ, u), ∂xᶜᶜᶜ(i,   j, k, grid, u), ∂xᶜᶜᶜ(i+1, j, k, grid, u))
    @inbounds s₂ = (∂xᶜᶜᶜ(i,   j, k, grid, ∂xᶜᶜᶜ, u), ∂xᶜᶜᶜ(i+1, j, k, grid, u), ∂xᶜᶜᶜ(i+2, j, k, grid, u))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_xᵃᶠᵃ(i, j, k, grid, v)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (v[i-2, j, k], v[i-1, j, k], v[i,   j, k])
    @inbounds ψ₁ = (v[i-1, j, k], v[i,   j, k], v[i+1, j, k])
    @inbounds ψ₂ = (v[i,   j, k], v[i+1, j, k], v[i+2, j, k])

    @inbounds s₀ = (∂yᶜᶜᶜ(i-2, j, k, grid, v), ∂yᶜᶜᶜ(i-1, j, k, grid, ∂yᶜᶜᶜ, v), ∂yᶜᶜᶜ(i,   j, k, grid, v))
    @inbounds s₁ = (∂yᶜᶜᶜ(i-1, j, k, grid, v), ∂yᶜᶜᶜ(i,   j, k, grid, ∂yᶜᶜᶜ, v), ∂yᶜᶜᶜ(i+1, j, k, grid, v))
    @inbounds s₂ = (∂yᶜᶜᶜ(i,   j, k, grid, v), ∂yᶜᶜᶜ(i+1, j, k, grid, ∂yᶜᶜᶜ, v), ∂yᶜᶜᶜ(i+2, j, k, grid, v))
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_yᵃᶠᵃ(i, j, k, grid, v)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (v[i, j-2, k], v[i, j-1, k], v[i, j,   k])
    @inbounds ψ₁ = (v[i, j-1, k], v[i, j,   k], v[i, j+1, k])
    @inbounds ψ₂ = (v[i, j,   k], v[i, j+1, k], v[i, j+2, k])

    @inbounds s₀ = (∂yᶜᶜᶜ(i, j-2, k, grid,  v), ∂yᶜᶜᶜ(i, j-1, k, grid,  v), ∂yᶜᶜᶜ(i, j,   k, grid,  v))
    @inbounds s₁ = (∂yᶜᶜᶜ(i, j-1, k, grid,  v), ∂yᶜᶜᶜ(i, j,   k, grid,  v), ∂yᶜᶜᶜ(i, j+1, k, grid,  v))
    @inbounds s₂ = (∂yᶜᶜᶜ(i, j,   k, grid,  v), ∂yᶜᶜᶜ(i, j+1, k, grid,  v), ∂yᶜᶜᶜ(i, j+2, k, grid,  v))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_yᶠᵃᵃ(i, j, k, grid, u)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (u[i, j-2, k], u[i, j-1, k], u[i, j,   k])
    @inbounds ψ₁ = (u[i, j-1, k], u[i, j,   k], u[i, j+1, k])
    @inbounds ψ₂ = (u[i, j,   k], u[i, j+1, k], u[i, j+2, k])

    @inbounds s₀ = (∂xᶜᶜᶜ(i, j-2, k, grid, u), ∂xᶜᶜᶜ(i, j-1, k, grid, u), ∂xᶜᶜᶜ(i, j,   k, grid, u))
    @inbounds s₁ = (∂xᶜᶜᶜ(i, j-1, k, grid, u), ∂xᶜᶜᶜ(i, j,   k, grid, u), ∂xᶜᶜᶜ(i, j+1, k, grid, u))
    @inbounds s₂ = (∂xᶜᶜᶜ(i, j,   k, grid, u), ∂xᶜᶜᶜ(i, j+1, k, grid, u), ∂xᶜᶜᶜ(i, j+2, k, grid, u))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end


@inline function center_interpolate_xᶠᵃᵃ(i, j, k, grid, f, args...)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (f(i-2, j, k, grid, args...), f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...))
    @inbounds ψ₁ = (f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...))
    @inbounds ψ₂ = (f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...), f(i+2, j, k, grid, args...))
    
    @inbounds s₀ = (∂xᶜᶜᶜ(i-2, j, k, grid, f, args...), ∂xᶜᶜᶜ(i-1, j, k, grid,  f, args...), ∂xᶜᶜᶜ(i,   j, k, grid,  f, args...))
    @inbounds s₁ = (∂xᶜᶜᶜ(i-1, j, k, grid, f, args...), ∂xᶜᶜᶜ(i,   j, k, grid,  f, args...), ∂xᶜᶜᶜ(i+1, j, k, grid,  f, args...))
    @inbounds s₂ = (∂xᶜᶜᶜ(i,   j, k, grid, f, args...), ∂xᶜᶜᶜ(i+1, j, k, grid,  f, args...), ∂xᶜᶜᶜ(i+2, j, k, grid,  f, args...))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_xᵃᶠᵃ(i, j, k, grid, f, args...)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (f(i-2, j, k, grid, args...), f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...))
    @inbounds ψ₁ = (f(i-1, j, k, grid, args...), f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...))
    @inbounds ψ₂ = (f(i,   j, k, grid, args...), f(i+1, j, k, grid, args...), f(i+2, j, k, grid, args...))

    @inbounds s₀ = (∂yᶜᶜᶜ(i-2, j, k, grid, f, args...), ∂yᶜᶜᶜ(i-1, j, k, grid, f, args...), ∂yᶜᶜᶜ(i,   j, k, grid, f, args...))
    @inbounds s₁ = (∂yᶜᶜᶜ(i-1, j, k, grid, f, args...), ∂yᶜᶜᶜ(i,   j, k, grid, f, args...), ∂yᶜᶜᶜ(i+1, j, k, grid, f, args...))
    @inbounds s₂ = (∂yᶜᶜᶜ(i,   j, k, grid, f, args...), ∂yᶜᶜᶜ(i+1, j, k, grid, f, args...), ∂yᶜᶜᶜ(i+2, j, k, grid, f, args...))
    
    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_yᵃᶠᵃ(i, j, k, grid, f, args...)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (f(i, j-2, k, grid, args...), f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...))
    @inbounds ψ₁ = (f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...))
    @inbounds ψ₂ = (f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...), f(i, j+2, k, grid, args...))

    @inbounds s₀ = (∂yᶜᶜᶜ(i, j-2, k, grid,  f, args...), ∂yᶜᶜᶜ(i, j-1, k, grid,  f, args...), ∂yᶜᶜᶜ(i, j,   k, grid,  f, args...))
    @inbounds s₁ = (∂yᶜᶜᶜ(i, j-1, k, grid,  f, args...), ∂yᶜᶜᶜ(i, j,   k, grid,  f, args...), ∂yᶜᶜᶜ(i, j+1, k, grid,  f, args...))
    @inbounds s₂ = (∂yᶜᶜᶜ(i, j,   k, grid,  f, args...), ∂yᶜᶜᶜ(i, j+1, k, grid,  f, args...), ∂yᶜᶜᶜ(i, j+2, k, grid,  f, args...))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

@inline function center_interpolate_yᶠᵃᵃ(i, j, k, grid, f, args...)
    
    FT = eltype(grid)

    @inbounds ψ₀ = (f(i, j-2, k, grid, args...), f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...))
    @inbounds ψ₁ = (f(i, j-1, k, grid, args...), f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...))
    @inbounds ψ₂ = (f(i, j,   k, grid, args...), f(i, j+1, k, grid, args...), f(i, j+2, k, grid, args...))

    @inbounds s₀ = (∂xᶜᶜᶜ(i, j-2, k, grid, f, args...), ∂xᶜᶜᶜ(i, j-1, k, grid, f, args...), ∂xᶜᶜᶜ(i, j,   k, grid, f, args...))
    @inbounds s₁ = (∂xᶜᶜᶜ(i, j-1, k, grid, f, args...), ∂xᶜᶜᶜ(i, j,   k, grid, f, args...), ∂xᶜᶜᶜ(i, j+1, k, grid, f, args...))
    @inbounds s₂ = (∂xᶜᶜᶜ(i, j,   k, grid, f, args...), ∂xᶜᶜᶜ(i, j+1, k, grid, f, args...), ∂xᶜᶜᶜ(i, j+2, k, grid, f, args...))

    return weno_reconstruction(FT, ψ₀, ψ₁, ψ₂, s₀, s₁, s₂)
end

using Oceananigans.Grids: inactive_node
using Oceananigans.ImmersedBoundaries: calc_inactive_stencil

const f = Face()
const c = Center()

@eval begin
    @inline near_boundary_center_xᶠᵃᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(4, :left, :x, :ᶜ, xside = :ᶜ)...))
    @inline near_boundary_center_yᵃᶠᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(4, :left, :y, :ᶜ, yside = :ᶜ)...))
    @inline near_boundary_center_xᵃᶠᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(3, :left, :x, :ᶜ, xside = :ᶜ)...), $(calc_inactive_stencil(3, :left, :x, :ᶜ, xside = :ᶜ, yshift = 1)...))
    @inline near_boundary_center_yᶠᵃᵃ(i, j, k, ibg) = @inbounds (|)($(calc_inactive_stencil(3, :left, :y, :ᶜ, yside = :ᶜ)...), $(calc_inactive_stencil(3, :left, :y, :ᶜ, yside = :ᶜ, xshift = 1)...))
end

@inline _center_interpolate_yᵃᶠᵃ(i, j, k, grid, v) = 
        ifelse(near_boundary_center_yᵃᶠᵃ(i, j, k, grid),  v[i, j, k], center_interpolate_yᵃᶠᵃ(i, j, k, grid, v))

@inline _center_interpolate_xᶠᵃᵃ(i, j, k, grid, u) = 
        ifelse(near_boundary_center_xᶠᵃᵃ(i, j, k, grid),  u[i, j, k], center_interpolate_xᶠᵃᵃ(i, j, k, grid, u))

@inline _center_interpolate_yᶠᵃᵃ(i, j, k, grid, v) = 
        ifelse(near_boundary_center_yᶠᵃᵃ(i, j, k, grid),  v[i, j, k], center_interpolate_yᶠᵃᵃ(i, j, k, grid, v))

@inline _center_interpolate_xᵃᶠᵃ(i, j, k, grid, u) = 
        ifelse(near_boundary_center_xᵃᶠᵃ(i, j, k, grid),  u[i, j, k], center_interpolate_xᵃᶠᵃ(i, j, k, grid, u))

@inline _center_interpolate_yᵃᶠᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_yᵃᶠᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_yᵃᶠᵃ(i, j, k, grid, f, args...))

@inline _center_interpolate_xᶠᵃᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_xᶠᵃᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_xᶠᵃᵃ(i, j, k, grid, f, args...))

@inline _center_interpolate_yᶠᵃᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_yᶠᵃᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_yᶠᵃᵃ(i, j, k, grid, f, args...))

@inline _center_interpolate_xᵃᶠᵃ(i, j, k, grid, f::Function, args...) = 
        ifelse(near_boundary_center_xᵃᶠᵃ(i, j, k, grid),  f(i, j, k, grid, args...), center_interpolate_xᵃᶠᵃ(i, j, k, grid, f, args...))

@inline flux_div_xyᶜᶜᶠ(i, j, k, grid, Qu, Qv) = δxᶜᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᶠᵃᵃ, Qu) + δyᵃᶜᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᵃᶠᵃ, Qv)

@inline div_xyᶜᶜᶜ(i, j, k, grid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᶠᵃᵃ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᵃᶠᵃ, v))

@inline div_xyᶜᶜᶜ(i, j, k, grid::ImmersedBoundaryGrid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᶠᵃᵃ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᵃᶠᵃ, v))

@inline div_xyᶜᶜᶠ(i, j, k, grid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶠ, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᶠᵃᵃ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶠ, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᵃᶠᵃ, v))

@inline div_xyᶜᶜᶠ(i, j, k, grid::ImmersedBoundaryGrid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶠ, _center_interpolate_xᶠᵃᵃ, _center_interpolate_yᶠᵃᵃ, u) +
                                δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶠ, _center_interpolate_yᵃᶠᵃ, _center_interpolate_xᵃᶠᵃ, v))
                            
@inline function advective_tracer_flux_x(i, j, k, grid, scheme::UpwindScheme, U, c) 

    ũ  =       _center_interpolate_yᶠᵃᵃ(i, j, k, grid, _center_interpolate_xᶠᵃᵃ, U)
    cᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, c)
    cᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, c)

    return Axᶠᶜᶜ(i, j, k, grid) * upwind_biased_product(ũ, cᴸ, cᴿ)
end

@inline function advective_tracer_flux_y(i, j, k, grid, scheme::UpwindScheme, V, c)

    ṽ  =       _center_interpolate_xᵃᶠᵃ(i, j, k, grid, _center_interpolate_yᵃᶠᵃ, V)
    cᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, c)
    cᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, c)

    return Ayᶜᶠᶜ(i, j, k, grid) * upwind_biased_product(ṽ, cᴸ, cᴿ)
end