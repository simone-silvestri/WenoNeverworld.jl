import Oceananigans.TurbulenceClosures:
        compute_diffusivities!,
        DiffusivityFields,
        viscosity, 
        diffusivity,
        diffusive_flux_x,
        diffusive_flux_y, 
        diffusive_flux_z,
        viscous_flux_ux,
        viscous_flux_uy,
        viscous_flux_uz,
        viscous_flux_vx,
        viscous_flux_vy,
        viscous_flux_vz,
        viscous_flux_wx,
        viscous_flux_wy,
        viscous_flux_wz

using Oceananigans.BuoyancyModels: ∂x_b, ∂y_b, ∂z_b 

"""
    struct EnergyBackScattering{FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 3}

Energy backscattering turbulence closure model.
This struct represents a turbulence closure model based on the energy backscattering principle. 
It is a parameterization of the turbulent momentum flux in a fluid flow. 
The model is implemented as a struct with a type parameter `FT` representing the floating-point type used for calculations.

# Arguments
- `ν::FT`: The kinematic anti-viscosity of the fluid.

reference:
Zanna, L., Bolton, T. (2020). 
Data-driven equation discovery of ocean mesoscale closures. 
Geophysical Research Letters, 47, e2020GL088376. https://doi.org/10.1029/2020GL088376
"""
struct EnergyBackScattering{FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 3}
    ν :: FT
end

EnergyBackScattering(FT::DataType = Float64; ν=FT(-4.87e7)) = EnergyBackScattering(ν) 

const MBS = EnergyBackScattering

@inline D̃ᶜᶜᶜ(i, j, k, grid, u, v) = 1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᶜᶜ(i, j, k, grid, Ax_qᶜᶜᶜ, u) -  
                                                               δyᶜᶜᶜ(i, j, k, grid, Ay_qᶜᶜᶜ, v))
                                                              
@inline Dᶠᶠᶜ(i, j, k, grid, u, v) = 1 / Vᶠᶠᶜ(i, j, k, grid) * (δyᶠᶠᶜ(i, j, k, grid, Ay_qᶠᶠᶜ, u) +  
                                                               δxᶠᶠᶜ(i, j, k, grid, Ax_qᶠᶠᶜ, v))
                                                              
#####
##### Abstract Smagorinsky functionality
#####

@inline ν(closure::MBS) = closure.ν

# Vertical viscous fluxes for isotropic diffusivities
@inline viscous_flux_uz(i, j, k, grid, clo::MBS, K, clk, fields, b) = zero(grid)
@inline viscous_flux_vz(i, j, k, grid, clo::MBS, K, clk, fields, b) = zero(grid)
@inline viscous_flux_wz(i, j, k, grid, clo::MBS, K, clk, fields, b) = zero(grid)
@inline viscous_flux_wx(i, j, k, grid, clo::MBS, K, clk, fields, b) = zero(grid)
@inline viscous_flux_wy(i, j, k, grid, clo::MBS, K, clk, fields, b) = zero(grid)

@inline ζ²_ζDᶠᶠᶜ(i, j, k, grid, u, v) = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) * (ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) - Dᶠᶠᶜ(i, j, k, grid, u, v))
@inline    ζD̃ᶠᶠᶜ(i, j, k, grid, u, v) = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) * ℑxyᶠᶠᵃ(i, j, k, grid, D̃ᶜᶜᶜ, u, v)

@inline viscous_flux_ux(i, j, k, grid, clo::MBS, K, clk, fields, b) = - ν(clo) * ℑxyᶜᶜᵃ(i, j, k, grid, ζ²_ζDᶠᶠᶜ, fields.u, fields.v)
@inline viscous_flux_vx(i, j, k, grid, clo::MBS, K, clk, fields, b) = - ν(clo) * ζD̃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
@inline viscous_flux_uy(i, j, k, grid, clo::MBS, K, clk, fields, b) = - ν(clo) * ζD̃ᶠᶠᶜ(i, j, k, grid, fields.u, fields.v)
@inline viscous_flux_vy(i, j, k, grid, clo::MBS, K, clk, fields, b) = - ν(clo) * ℑxyᶜᶜᵃ(i, j, k, grid, ζ²_ζDᶠᶠᶜ, fields.u, fields.v)

@inline diffusive_flux_x(i, j, k, grid, closure::MBS, K, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::MBS, K, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::MBS, K, ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index = zero(grid)
