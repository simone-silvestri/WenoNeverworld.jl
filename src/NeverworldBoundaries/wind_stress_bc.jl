struct WindStressBoundaryCondition{F, T, S} <: Function
    φs :: F
    τs :: T
    stress :: S
end

default_φs = (-70, -45, -15, 0, 15, 45, 70)
default_τs = (0.0, 0.2, -0.1, -0.02, -0.1, 0.1, 0.0)
    
"""
    WindStressBoundaryCondition(; φs = default_φs, τs = default_τs)

Wind stess boundary condition which implements a piecewise cubic interpolation
between points `φs` (`Tuple`) and `τs` (`Tuple`).
"""
WindStressBoundaryCondition(; φs = default_φs, τs = default_τs) =  WindStressBoundaryCondition(φs, τs, nothing)

(ws::WindStressBoundaryCondition)(i, j, grid, clock, fields) = ws.stress[j]

Adapt.adapt_structure(to, ws::WindStressBoundaryCondition) = WindStressBoundaryCondition(nothing, nothing, adapt(to, ws.stress))

@inline function regularize_boundary_condition(bc::WindStressBoundaryCondition, grid)

    Ny   = size(grid, 2)
    arch = architecture(grid)
    
    φ_grid = grid.φᵃᶜᵃ[1:Ny]

    stress = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)    
        φ_index = sum(φ .> bc.φs) + 1
        φ₁ = bc.φs[φ_index-1]
        φ₂ = bc.φs[φ_index]
        τ₁ = bc.τs[φ_index-1]
        τ₂ = bc.τs[φ_index]
        stress[j] = cubic_interpolate(φ, x₁ = φ₁, x₂ = φ₂, y₁ = τ₁, y₂ = τ₂) / 1000.0
    end

    return WindStressBoundaryCondition(bc.φs, bc.τs, on_architecture(arch, - stress))
end