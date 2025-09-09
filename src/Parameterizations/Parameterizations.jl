module Parameterizations

export QGLeith, EnergyBackScattering, OMp25Closure

using Oceananigans
using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: 
        AbstractTurbulenceClosure,
        AbstractScalarBiharmonicDiffusivity,
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
        isopycnal_rotation_tensor_xz_ccf,
        isopycnal_rotation_tensor_yz_ccf,
        isopycnal_rotation_tensor_zz_ccf

import Oceananigans.TurbulenceClosures:
        compute_diffusivities!,
        build_diffusivity_fields,
        viscosity, 
        diffusivity,
        diffusive_flux_x,
        diffusive_flux_y, 
        diffusive_flux_z

using Oceananigans.Utils: launch!
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Operators
using Oceananigans.BuoyancyFormulations: ∂x_b, ∂y_b, ∂z_b 

using Oceananigans.Operators: ℑxyzᶜᶜᶠ, ℑyzᵃᶜᶠ, ℑxzᶜᵃᶠ, Δxᶜᶜᶜ, Δyᶜᶜᶜ
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz

"The averaged filter width"
@inline Δ̃ᶜᶜᶜ(i, j, k, grid) = sqrt((Δxᶜᶜᶜ(i, j, k, grid)^2 + Δyᶜᶜᶜ(i, j, k, grid)^2)/2)

"Return the filter width for an Horizontal closure on a general grid."
@inline Δ²ᶜᶜᶜ(i, j, k, grid) =  2 * (1 / (1 / Δxᶜᶜᶜ(i, j, k, grid)^2 + 1 / Δyᶜᶜᶜ(i, j, k, grid)^2))

include("quasi_geostrophic_leith.jl")
include("energy_backscattering.jl")
include("smagorinsky_closure.jl")

end