module UpwindVectorInvariantSchemes

export UpwindVectorInvariant, GlobalVectorInvariant, UpwindKineticScheme, UpwindVerticalScheme
export WENO, VelocityStencil, VorticityStencil

using Oceananigans
using WenoNeverworld
using Oceananigans.Advection
using Adapt

using Oceananigans.Advection: AbstractAdvectionScheme
import Oceananigans.Advection: WENO, VorticityStencil, VelocityStencil

import Oceananigans.Advection: vertical_advection_U, vertical_advection_V
import Oceananigans.Advection: bernoulli_head_U, bernoulli_head_V
import Oceananigans.Advection: U_dot_∇u, U_dot_∇v
import Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

using Oceananigans.Advection: _advective_momentum_flux_Wu, 
                              _advective_momentum_flux_Wv,
                              _advective_momentum_flux_Uu,
                              _advective_momentum_flux_Vv

using Oceananigans.Advection:  _left_biased_interpolate_xᶠᵃᵃ,
                               _left_biased_interpolate_yᵃᶠᵃ,
                              _right_biased_interpolate_xᶠᵃᵃ,
                              _right_biased_interpolate_yᵃᶠᵃ,
                               _left_biased_interpolate_xᶜᵃᵃ,
                               _left_biased_interpolate_yᵃᶜᵃ,
                              _right_biased_interpolate_xᶜᵃᵃ,
                              _right_biased_interpolate_yᵃᶜᵃ

using Oceananigans.Advection: upwind_biased_product
using Oceananigans.Advection: compute_reconstruction_coefficients, SmoothnessStencil

using Oceananigans.Operators

import Oceananigans.Models.HydrostaticFreeSurfaceModels: validate_momentum_advection

@inline smoothness_stencil(scheme) = VorticityStencil

include("extend_weno_properties.jl")
include("multi_dimensional_reconstrution.jl")
include("default_vector_invariant.jl")
include("vertical_advection_treatment.jl")
include("global_upwind_vector_invariant.jl")

end 