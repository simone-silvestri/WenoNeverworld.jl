using Oceananigans.TurbulenceClosures
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using CUDA: @allowscalar 

@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Δᶜᶜᶜ(i, j, k, grid)    = min(Δxᶜᶜᶜ(i, j, k, grid), Δyᶜᶜᶜ(i, j, k, grid))

@inline function νhb_smagorinski_final(i, j, k, grid, clock, fields, C₄)
    δ₁ = Dₛ(i, j, k, grid, fields.u, fields.v)
    δ₂ = ℑxyᶜᶜᵃ(i, j, k, grid, Dₜ, fields.u, fields.v)
    return Δᶜᶜᶜ(i, j, k, grid)^4 * C₄ * sqrt(δ₁^2 + δ₂^2)
end

function smagorinsky_viscosity(formulation, grid; Cₛₘ = 4.0)

    dx_min = min_Δx(grid.underlying_grid)
    dy_min = min_Δy(grid.underlying_grid)
    dx_max = @allowscalar grid.Δxᶠᶜᵃ[Int(grid.Ny / 2)]
    dy_max = grid.Δyᶠᶜᵃ
    timescale_max = 100days
    timescale_min = 0.2days

    @show C₄    = (Cₛₘ / π)^2 / 8
    @show min_ν = (1 / (1 / dx_min^2 + 1 / dy_min^2))^2 / timescale_max
    @show max_ν = (1 / (1 / dx_max^2 + 1 / dy_max^2))^2 / timescale_min

    loc = (Center, Center, Center)

    return ScalarBiharmonicDiffusivity(formulation;
                                       ν=νhb_smagorinski_final, discrete_form=true, loc,
                                       parameters = C₄)
end

@inline function νhb_leith_final(i, j, k, grid, clock, fields, p) 
    ∂ζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)^2 + ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)^2
    ∂δ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)^2 + ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)^2

    return max(Δᶜᶜᶜ(i, j, k, grid)^5 * p.C₄ * sqrt(∂ζ + ∂δ), p.min_ν)
end

function leith_viscosity(formulation, grid; Cₗₜ = 1.0)

    dx_min    = min_Δx(grid.underlying_grid)
    dy_min    = min_Δy(grid.underlying_grid)
    timescale = 10days

    @show C₄    = (Cₗₜ / π)^3 / 8
    @show min_ν = (1 / (1 / dx_min^2 + 1 / dy_min^2))^2 / timescale

    loc = (Center, Center, Center)

    return ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_leith_final, discrete_form=true, loc, 
                                       parameters = (C₄ = C₄, min_ν = min_ν))
end

using Oceananigans.Operators: Δx, Δy

@inline νhb_geometric(i, j, k, grid, lx, ly, lz, clock, fields, λ) =
                      (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))^2 / λ

geometric_viscosity(formulation, timescale) = ScalarBiharmonicDiffusivity(formulation, ν=νhb_geometric, discrete_form=true, parameters = timescale) 
