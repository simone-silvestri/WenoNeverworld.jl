using WenoNeverworld.Auxiliaries: parabolic_scaling

struct BuoyancyRelaxationBoundaryCondition{T, S, F} <: Function
    ΔB::T
    λ::S
    func::F
end

"""
    BuoyancyRelaxationBoundaryCondition(func = (y, t) -> parabolic_scaling(y); ΔB = ΔB, λ = 7days)

Buoyancy relaxation profile which implements a latitude-time dependent boundary condition following: 

`b = Δz_surface / λ * (b_surf - ΔB * func(φ, t))`

Arguments:
==========

- func: function which takes the latitude φ and time t and returns a scalar

Keyword arguments:
==================

- ΔB: buoyancy difference between the equator and the poles, default: 6.0e-2
- λ: restoring time-scale, default: 7days
"""
BuoyancyRelaxationBoundaryCondition(func = (y, t) -> parabolic_scaling(y); ΔB = Constants.ΔB, λ = 7days) = BuoyancyRelaxationBoundaryCondition(ΔB, λ, func)

function (b::BuoyancyRelaxationBoundaryCondition)(i, j, grid, clock, fields)
    φ  = φnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    Δz = Δzᶜᶜᶜ(i, j, grid.Nz, grid)
    b_surf = fields.b[i, j, grid.Nz]
    return Δz / b.λ * (b_surf - b.ΔB * b.func(φ, clock.time))
end

Adapt.adapt_structure(to, b::BuoyancyRelaxationBoundaryCondition) = BuoyancyRelaxationBoundaryCondition(b.ΔB, b.λ, b.func)
