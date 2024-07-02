module Parameterizations

export QGLeith, 
       EnergyBackScattering, 
       NNbackscatteringClosure,
       XinKaiVerticalDiffusivity

using Oceananigans
using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: 
        AbstractTurbulenceClosure,
        HorizontalFormulation,
        HorizontalDivergenceFormulation, 
        HorizontalDivergenceScalarBiharmonicDiffusivity

using Oceananigans.TurbulenceClosures:
        tapering_factorᶠᶜᶜ,
        tapering_factorᶜᶠᶜ,
        tapering_factorᶜᶜᶠ,
        tapering_factor,
        SmallSlopeIsopycnalTensor,
        AbstractScalarDiffusivity,
        VerticallyImplicitTimeDiscretization,
        ExplicitTimeDiscretization,
        FluxTapering,
        getclosure,
        isopycnal_rotation_tensor_xz_ccf,
        isopycnal_rotation_tensor_yz_ccf,
        isopycnal_rotation_tensor_zz_ccf

import Oceananigans.TurbulenceClosures:
        compute_diffusivities!,
        DiffusivityFields,
        viscosity, 
        diffusivity,
        diffusive_flux_x,
        diffusive_flux_y, 
        diffusive_flux_z,
        top_buoyancy_flux

import Oceananigans.TurbulenceClosures: 
        ∂ⱼ_τ₁ⱼ, 
        ∂ⱼ_τ₂ⱼ, 
        ∂ⱼ_τ₃ⱼ,
        ∇_dot_qᶜ

using Oceananigans.Utils: launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Operators
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BuoyancyModels: ∂x_b, ∂y_b, ∂z_b 

using Oceananigans.Operators: ℑxyzᶜᶜᶠ, ℑyzᵃᶜᶠ, ℑxzᶜᵃᶠ, Δxᶜᶜᶜ, Δyᶜᶜᶜ

using Adapt

"Return the filter width for an Horizontal closure on a general grid."
@inline Δ²ᶜᶜᶜ(i, j, k, grid) =  2 * (1 / (1 / Δxᶜᶜᶜ(i, j, k, grid)^2 + 1 / Δyᶜᶜᶜ(i, j, k, grid)^2))

include("quasi_geostrophic_leith.jl")
include("energy_backscattering.jl")
include("xin_kai_vertical_diffusivity.jl")
include("guillaumin_zanna_parameterization.jl")

end