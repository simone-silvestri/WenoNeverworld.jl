using Oceananigans.Advection: AbstractUpwindBiasedAdvectionScheme,
                                      _symmetric_interpolate_xᶠᵃᵃ,
                                      _symmetric_interpolate_yᵃᶠᵃ,
                                      _symmetric_interpolate_zᵃᵃᶠ,
                                    _left_biased_interpolate_xᶠᵃᵃ,
                                    _left_biased_interpolate_yᵃᶠᵃ,
                                    _left_biased_interpolate_zᵃᵃᶠ,
                                   _right_biased_interpolate_xᶠᵃᵃ,
                                   _right_biased_interpolate_yᵃᶠᵃ,
                                   _right_biased_interpolate_zᵃᵃᶠ,
                                            upwind_biased_product,
                                          advective_tracer_flux_x,
                                          advective_tracer_flux_y,
                                          advective_tracer_flux_z
                                            

import Oceananigans.Advection: div_Uc

using Oceananigans.TurbulenceClosures: SmallSlopeIsopycnalTensor, FluxTapering
using Oceananigans.Operators
using Adapt

abstract type AbstractIsopycnallyRotatedUpwindBiasedAdvection{N, FT} <: AbstractUpwindBiasedAdvectionScheme{N, FT} end

struct IsopycnallyRotatedUpwindScheme{N, FT, U, C, I, S} <: AbstractIsopycnallyRotatedUpwindBiasedAdvection{N, FT}
    upwind_scheme    :: U
    centered_scheme  :: C
    isopycnal_tensor :: I
    slope_limiter    :: S

    IsopycnallyRotatedUpwindScheme{N, FT}(u::U, c::C, i::I, s::S) where {N, FT, U, C, I, S} = new{N, FT, U, C, I, S}(u, c, i, s)
end

function IsopycnallyRotatedUpwindScheme(upwind_scheme::AbstractUpwindBiasedAdvectionScheme{N, FT}, centered_scheme; 
                                        isopycnal_tensor = SmallSlopeIsopycnalTensor(),
                                        slope_limiter = FluxTapering(1e-2)) where {N, FT}

    return IsopycnallyRotatedUpwindScheme{N, FT}(upwind_scheme, centered_scheme, isopycnal_tensor, slope_limiter)
end

Adapt.adapt_structure(to, scheme::IsopycnallyRotatedUpwindScheme{N, FT}) where {N, FT} =
    IsopycnallyRotatedUpwindScheme{N, FT}(Adapt.adapt(to, scheme.upwind_scheme),    Adapt.adapt(to, scheme.centered_scheme),
                                          Adapt.adapt(to, scheme.isopycnal_tensor), Adapt.adapt(to, scheme.slope_limiter))


import Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_free_surface_tracer_tendency
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ∇_dot_qᶜ, immersed_∇_dot_qᶜ, hydrostatic_fields

@inline function hydrostatic_free_surface_tracer_tendency(i, j, k, grid,
                                                          val_tracer_index::Val{tracer_index},
                                                          advection::IsopycnallyRotatedUpwindScheme,
                                                          closure,
                                                          c_immersed_bc,
                                                          buoyancy,
                                                          velocities,
                                                          free_surface,
                                                          tracers,
                                                          top_tracer_bcs,
                                                          diffusivities,
                                                          auxiliary_fields,
                                                          forcing,
                                                          clock) where tracer_index

    @inbounds c = tracers[tracer_index]
    @inbounds b = tracers.b
    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    
    return ( - div_Uc(i, j, k, grid, advection, velocities, c, b)
             - ∇_dot_qᶜ(i, j, k, grid, closure, diffusivities, val_tracer_index, c, clock, model_fields, buoyancy)
             - immersed_∇_dot_qᶜ(i, j, k, grid, c, c_immersed_bc, closure, diffusivities, val_tracer_index, clock, model_fields)
             + forcing(i, j, k, grid, clock, model_fields))
end

@inline function div_Uc(i, j, k, grid, advection::AbstractIsopycnallyRotatedUpwindBiasedAdvection, U, c, b)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_x, advection, U, c, b) +
                                    δyᵃᶜᵃ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_y, advection, U, c, b) +
                                    δzᵃᵃᶜ(i, j, k, grid, _isopycnally_rotated_advective_tracer_flux_z, advection, U, c, b))
end
    
using Oceananigans.ImmersedBoundaries: conditional_flux_fcc, conditional_flux_cfc, conditional_flux_ccf, GFIBG
using Oceananigans.ImmersedBoundaries: near_x_immersed_boundary_symmetricᶜ, 
                                       near_y_immersed_boundary_symmetricᶜ, 
                                       near_z_immersed_boundary_symmetricᶜ,
                                       near_x_immersed_boundary_symmetricᶠ, 
                                       near_y_immersed_boundary_symmetricᶠ, 
                                       near_z_immersed_boundary_symmetricᶠ

@inline near_any_boundary(i, j, k, ibg, advection::IsopycnallyRotatedUpwindScheme) = 
                        near_x_immersed_boundary_symmetricᶜ(i, j, k, ibg, advection.centered_scheme) |
                        near_y_immersed_boundary_symmetricᶜ(i, j, k, ibg, advection.centered_scheme) |
                        near_z_immersed_boundary_symmetricᶜ(i, j, k, ibg, advection.centered_scheme) |
                        near_x_immersed_boundary_symmetricᶠ(i, j, k, ibg, advection.centered_scheme) |
                        near_y_immersed_boundary_symmetricᶠ(i, j, k, ibg, advection.centered_scheme) |
                        near_z_immersed_boundary_symmetricᶠ(i, j, k, ibg, advection.centered_scheme) 

@inline _isopycnally_rotated_advective_tracer_flux_x(args...) = isopycnally_rotated_advective_tracer_flux_x(args...)
@inline _isopycnally_rotated_advective_tracer_flux_y(args...) = isopycnally_rotated_advective_tracer_flux_y(args...)
@inline _isopycnally_rotated_advective_tracer_flux_z(args...) = isopycnally_rotated_advective_tracer_flux_z(args...)

@inline _isopycnally_rotated_advective_tracer_flux_x(i, j, k, ibg::GFIBG, adv, U, c, b) = ifelse(near_any_boundary(i, j, k, ibg, adv), _advective_tracer_flux_x(i, j, k, ibg, adv.upwind_scheme, U.u, c), isopycnally_rotated_advective_tracer_flux_x(i, j, k, ibg.underlying_grid, adv, U, c, b))
@inline _isopycnally_rotated_advective_tracer_flux_y(i, j, k, ibg::GFIBG, adv, U, c, b) = ifelse(near_any_boundary(i, j, k, ibg, adv), _advective_tracer_flux_y(i, j, k, ibg, adv.upwind_scheme, U.v, c), isopycnally_rotated_advective_tracer_flux_y(i, j, k, ibg.underlying_grid, adv, U, c, b))
@inline _isopycnally_rotated_advective_tracer_flux_z(i, j, k, ibg::GFIBG, adv, U, c, b) = ifelse(near_any_boundary(i, j, k, ibg, adv), _advective_tracer_flux_z(i, j, k, ibg, adv.upwind_scheme, U.w, c), isopycnally_rotated_advective_tracer_flux_z(i, j, k, ibg.underlying_grid, adv, U, c, b))

@inline function advective_tracer_fluxes_x(i, j, k, grid, scheme::IsopycnallyRotatedUpwindScheme, U, c) 

    @inbounds ũ = U[i, j, k]

    cˢ =    _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Axᶠᶜᶜ(i, j, k, grid) * upwind_biased_product(ũ, cᴸ, cᴿ)
    advective_flux_center = Axᶠᶜᶜ(i, j, k, grid) * ũ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function advective_tracer_fluxes_y(i, j, k, grid, scheme::IsopycnallyRotatedUpwindScheme, V, c) 

    @inbounds ṽ = V[i, j, k]

    cˢ =    _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Ayᶜᶠᶜ(i, j, k, grid) * upwind_biased_product(ṽ, cᴸ, cᴿ)
    advective_flux_center = Ayᶜᶠᶜ(i, j, k, grid) * ṽ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function advective_tracer_fluxes_z(i, j, k, grid, scheme::IsopycnallyRotatedUpwindScheme, W, c) 

    @inbounds w̃ = W[i, j, k]

    cˢ =    _symmetric_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.centered_scheme, c)
    cᴸ =  _left_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.upwind_scheme,   c)
    cᴿ = _right_biased_interpolate_zᵃᵃᶠ(i, j, k, grid, scheme.upwind_scheme,   c)

    advective_flux_upwind = Azᶜᶜᶠ(i, j, k, grid) * upwind_biased_product(w̃, cᴸ, cᴿ)
    advective_flux_center = Azᶜᶜᶠ(i, j, k, grid) * w̃ * cˢ

    return advective_flux_upwind, advective_flux_center
end

@inline function isopycnally_rotated_advective_tracer_flux_x(i, j, k, grid, advection, U, c, b)
    
    R₁₁ = one(grid)
    R₁₃ = isopycnal_rotation_tensor_xz_fcc(i, j, k, grid, b, advection.isopycnal_tensor)

    Aᵘˣ, Aᶜˣ = advective_tracer_fluxes_x(i, j, k, grid, advection, U.u, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_fluxes_z(i, j, k, grid, advection, U.w, c)

    Dᵘˣ = Aᵘˣ - Aᶜˣ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶠᶜᶜ(i, j, k, grid, advection, c)

    return (R₁₁ * Dᵘˣ + R₁₃ * Dᵘᶻ + Aᶜˣ) * ϵ + Aᵘˣ *  (1 - ϵ)
end

@inline function isopycnally_rotated_advective_tracer_flux_y(i, j, k, grid, advection, U, c, b)
    
    R₂₂ = one(grid)
    R₂₃ = isopycnal_rotation_tensor_yz_cfc(i, j, k, grid, b, advection.isopycnal_tensor)

    Aᵘʸ, Aᶜʸ = advective_tracer_fluxes_y(i, j, k, grid, advection, U.v, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_fluxes_z(i, j, k, grid, advection, U.w, c)

    Dᵘʸ = Aᵘʸ - Aᶜʸ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶜᶠᶜ(i, j, k, grid, advection, b)

    return  (R₂₂ * Dᵘʸ + R₂₃ * Dᵘᶻ + Aᶜʸ) * ϵ + Aᵘʸ * (1 - ϵ)
end

@inline function isopycnally_rotated_advective_tracer_flux_z(i, j, k, grid, advection, U, c, b)

    R₃₁ = isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, b, advection.isopycnal_tensor)
    R₃₂ = isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, b, advection.isopycnal_tensor)
    R₃₃ = isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, b, advection.isopycnal_tensor)

    Aᵘˣ, Aᶜˣ = advective_tracer_fluxes_x(i, j, k, grid, advection, U.u, c)
    Aᵘʸ, Aᶜʸ = advective_tracer_fluxes_y(i, j, k, grid, advection, U.v, c)
    Aᵘᶻ, Aᶜᶻ = advective_tracer_fluxes_z(i, j, k, grid, advection, U.w, c)

    Dᵘˣ = Aᵘˣ - Aᶜˣ
    Dᵘʸ = Aᵘʸ - Aᶜʸ
    Dᵘᶻ = Aᵘᶻ - Aᶜᶻ

    ϵ = tapering_factorᶜᶜᶠ(i, j, k, grid, advection, b)

    return ϵ * (R₃₁ * Dᵘˣ +
                R₃₂ * Dᵘʸ +
                R₃₃ * Dᵘᶻ + Aᶜᶻ) + (1 - ϵ) * Aᵘᶻ
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

@inline function isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, b, slope_model::SmallSlopeIsopycnalTensor)

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
using Oceananigans.TurbulenceClosures: calc_tapering

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
