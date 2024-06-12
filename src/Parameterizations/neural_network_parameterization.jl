import Oceananigans.TurbulenceClosures: 
                        ∂ⱼ_τ₁ⱼ, 
                        ∂ⱼ_τ₂ⱼ, 
                        ∂ⱼ_τ₃ⱼ,
                        ∇_dot_qᶜ

const RequiredHalo = 5

# TD can be VerticallyImplicitTimeDiscretization, 
struct NNParameterization{M, C2, C3} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, RequiredHalo}
    matrix :: M
    param2 :: C2
    param3 :: C3
end

@inline function ∂ⱼ_τ₁ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy)
    @inbounds uᵢ   = model_fields.u[i,   j, k]
    @inbounds uᵢ₊₁ = model_fields.u[i+1, j, k]

    # compute your 
    return
end

@inline function ∂ⱼ_τ₂ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy)
    
end

@inline function ∂ⱼ_τ₃ⱼ(i, j, k, grid, closure, diffusivities, clock, model_fields, buoyancy)
    disc = time_discretization(closure)
    return 1 / Vᶜᶜᶠ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Ax_qᶠᶜᶠ, _viscous_flux_wx, disc, closure, args...) +
                                      δyᵃᶜᵃ(i, j, k, grid, Ay_qᶜᶠᶠ, _viscous_flux_wy, disc, closure, args...) +
                                      δzᵃᵃᶠ(i, j, k, grid, Az_qᶜᶜᶜ, _viscous_flux_wz, disc, closure, args...))
end

@inline ∇_dot_qᶜ(i, j, k, grid, closure::NNParameterization, args...) = zero(grid)

