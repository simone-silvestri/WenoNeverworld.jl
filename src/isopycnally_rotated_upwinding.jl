
coefficients_C6 = (2, -16, 74, 74, -16, 2.0) ./ 120

abstract type AbstractIsopycnallyRotatedUpwindBiasedAdvection{N, FT} <: AbstractUpwindBiasedAdvectionScheme{N, FT} end

struct IsopycnallyRotatedUpwindScheme{N, FT, U, C, I, S} <: AbstractIsopycnallyRotatedUpwindBiasedAdvection{N, FT}
    upwind_scheme    :: U
    centered_scheme  :: C
    isopycnal_tensor :: I
    slope_limiter    :: S
end

import Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y, _advective_tracer_flux_z

@inline function div_Uc(i, j, k, grid, advection::AbstractIsopycnallyRotatedUpwindBiasedAdvection, U, c)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_x, advection, U, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_y, advection, U, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_z, advection, U, c))
end
    
using Oceananigans.ImmersedBoundaries: conditional_flux_fcc, conditional_flux_cfc, conditional_flux_ccf, GFIBG

@inline _isopycnally_rotated_advective_tracer_flux_x(i, j, k, ibg::GFIBG, args...) = conditional_flux_fcc(i, j, k, ibg, zero(eltype(ibg)), isopycnally_rotated_advective_tracer_flux_x(i, j, k, ibg, args...))
@inline _isopycnally_rotated_advective_tracer_flux_y(i, j, k, ibg::GFIBG, args...) = conditional_flux_cfc(i, j, k, ibg, zero(eltype(ibg)), isopycnally_rotated_advective_tracer_flux_y(i, j, k, ibg, args...))
@inline _isopycnally_rotated_advective_tracer_flux_z(i, j, k, ibg::GFIBG, args...) = conditional_flux_ccf(i, j, k, ibg, zero(eltype(ibg)), isopycnally_rotated_advective_tracer_flux_z(i, j, k, ibg, args...))

@inline function advective_tracer_flux_x(i, j, k, grid, scheme::IsopycnallyRotatedUpwindBiasedAdvection, U, c) 

    @inbounds ũ = U[i, j, k]

    cˢ =    _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Axᶠᶜᶜ(i, j, k, grid) * upwind_biased_product(ũ, cᴸ, cᴿ)
    advective_flux_center = Axᶠᶜᶜ(i, j, k, grid) * ũ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function advective_tracer_flux_y(i, j, k, grid, scheme::IsopycnallyRotatedUpwindBiasedAdvection, V, c) 

    @inbounds ṽ = V[i, j, k]

    cˢ =    _symmetric_interpolate_yᶠᵃᵃ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_yᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_yᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Ayᶜᶠᶜ(i, j, k, grid) * upwind_biased_product(ṽ, cᴸ, cᴿ)
    advective_flux_center = Ayᶜᶠᶜ(i, j, k, grid) * ṽ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function advective_tracer_flux_z(i, j, k, grid, scheme::IsopycnallyRotatedUpwindBiasedAdvection, W, c) 

    @inbounds w̃ = W[i, j, k]

    cˢ =    _symmetric_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Azᶜᶜᶠ(i, j, k, grid) * upwind_biased_product(w̃, cᴸ, cᴿ)
    advective_flux_center = Azᶜᶜᶠ(i, j, k, grid) * w̃ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function isopycnally_rotated_advective_tracer_flux_x(i, j, k, grid, advection, U, c)
    
    R₁₁ = one(grid)
    R₁₃ = isopycnal_rotation_tensor_xz_fcc(i, j, k, grid, c, advection.isopycnal_tensor)

    Aᵘˣ, Aᶜˣ = advective_tracer_flux_x(i, j, k, grid, advection, U.u, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_flux_z(i, j, k, grid, advection, U.w, c)

    Dᵘˣ = Aᵘˣ - Aᶜˣ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶠᶜᶜ(i, j, k, grid, advection, c)

    return (R₁₁ * Dᵘˣ + R₁₃ * Dᵘᶻ + Aᶜˣ) * ϵ + Aᵘˣ *  (1 - ϵ)
end

@inline function isopycnally_rotated_advective_tracer_flux_y(i, j, k, grid, advection, U, c)
    
    R₂₂ = one(grid)
    R₂₃ = isopycnal_rotation_tensor_yz_cfc(i, j, k, grid, c, advection.isopycnal_tensor)

    Aᵘʸ, Aᶜʸ = advective_tracer_flux_y(i, j, k, grid, advection, U.v, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_flux_z(i, j, k, grid, advection, U.w, c)

    Dᵘʸ = Aᵘʸ - Aᶜʸ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶜᶠᶜ(i, j, k, grid, advection, c)

    return  (R₂₂ * Dᵘʸ + R₂₃ * Dᵘᶻ + Aᶜʸ) * ϵ + Aᵘʸ * (1 - ϵ)
end

@inline function isopycnally_rotated_advective_tracer_flux_z(i, j, k, grid, advection, U, c)

    R₃₁ = isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, c, advection.isopycnal_tensor)
    R₃₂ = isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, c, advection.isopycnal_tensor)
    R₃₃ = isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, c, advection.isopycnal_tensor)

    Aᵘˣ, Aᶜˣ = advective_tracer_flux_x(i, j, k, grid, advection, U.u, c)
    Aᵘʸ, Aᶜʸ = advective_tracer_flux_y(i, j, k, grid, advection, U.v, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_flux_z(i, j, k, grid, advection, U.w, c)

    Dᵘˣ = Aᵘˣ - Aᶜˣ
    Dᵘʸ = Aᵘʸ - Aᶜʸ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶜᶜᶠ(i, j, k, grid, advection, c)

    return ϵ * (R₃₁ * Dᵘˣ  +
                R₃₂ * Dᵘʸ  +
                R₃₃ * Dᵘᶻ  + Aᶜᶻ) + (1 - ϵ) * Aᵘᶻ
end

import Oceananigans.TurbulenceClosures: isopycnal_rotation_tensor_xz_fcc, 
                                        isopycnal_rotation_tensor_xz_ccf,
                                        isopycnal_rotation_tensor_yz_cfc,
                                        isopycnal_rotation_tensor_yz_ccf,
                                        isopycnal_rotation_tensor_zz_ccf
                                        
@inline function isopycnal_rotation_tensor_xz_fcc(i, j, k, grid, b, slope_model::SmallSlopeIsopycnalTensor)
    
    bx = ∂xᶠᶜᶜ(i, j, k, grid, b)
    bz = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
    bz = max(bz, slope_model.minimum_bz)

    slope_x = - bx / bz

    return ifelse(bz == 0, zero(grid), slope_x)
end

@inline function isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, b, slope_model::SmallSlopeIsopycnalTensor)

    bx = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    bz = max(bz, slope_model.minimum_bz)

    slope_x = - bx / bz

    return ifelse(bz == 0, zero(grid), slope_x)
end

@inline function isopycnal_rotation_tensor_yz_cfc(i, j, k, grid, b, slope_model::SmallSlopeIsopycnalTensor)

    by = ∂yᶜᶠᶜ(i, j, k, grid, b)
    bz = ℑyzᵃᶠᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
    bz = max(bz, slope_model.minimum_bz)

    slope_y = - by / bz

    return ifelse(bz == 0, zero(grid), slope_y)
end

@inline function isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, b, slope_model::SmallSlopeIsopycnalTensor)

    by = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    bz = max(bz, slope_model.minimum_bz)

    slope_y = - by / bz

    return ifelse(bz == 0, zero(grid), slope_y)
end

@inline function isopycnal_rotation_tensor_zz_ccf(i, j, k, grid::AbstractGrid, b, slope_model::SmallSlopeIsopycnalTensor)

    bx = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    by = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, b)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, b)
    bz = max(bz, slope_model.minimum_bz)

    slope_x = - bx / bz
    slope_y = - by / bz
    slope² = slope_x^2 + slope_y^2

    return ifelse(bz == 0, zero(grid), slope²)
end

import Oceananigans.TurbulenceClosures: tapering_factorᶠᶜᶜ, tapering_factorᶜᶠᶜ, tapering_factorᶜᶜᶠ


@inline function tapering_factorᶠᶜᶜ(i, j, k, grid, advection, c)
    
    by = ℑxyᶠᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, c)
    bz = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, c)
    bx = ∂xᶠᶜᶜ(i, j, k, grid, c)

    return calc_tapering(bx, by, bz, grid, advection.isopycnal_tensor, advection.slope_limiter)
end

@inline function tapering_factorᶜᶠᶜ(i, j, k, grid, advection, c)

    bx = ℑxyᶜᶠᵃ(i, j, k, grid, ∂xᶠᶜᶜ, c)
    bz = ℑyzᵃᶠᶜ(i, j, k, grid, ∂zᶜᶜᶠ, c)
    by = ∂yᶜᶠᶜ(i, j, k, grid, c)

    return calc_tapering(bx, by, bz, grid, advection.isopycnal_tensor, advection.slope_limiter)
end

@inline function tapering_factorᶜᶜᶠ(i, j, k, grid, advection, c)

    bx = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, c)
    by = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, c)
    bz = ∂zᶜᶜᶠ(i, j, k, grid, c)

    return calc_tapering(bx, by, bz, grid, advection.isopycnal_tensor, advection.slope_limiter)
end
