using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, AbstractScalarBiharmonicDiffusivity

import Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(ρ₀ * (1 - g * b)))

function HeightField(grid, loc = (Center, Center, Center))  

    zf = Field(loc, grid)

    for k in 1:size(zf, 3)
        interior(zf, :, :, k) .= znode(loc[3](), k, grid)
    end

    return zf
end

VerticalVorticityField(fields::Dict, i) = VerticalVorticityField((; u = fields[:u][i], v = fields[:v][i]))
KineticEnergyField(fields::Dict, i)     = KineticEnergyField((; u = fields[:u][i], v = fields[:v][i]))

function KineticEnergyField(velocities::NamedTuple)
    u = velocities.u
    v = velocities.v

    E_op = @at (Center, Center, Center) 0.5 * (u^2 + v^2)

    return compute!(Field(E_op))
end

function VerticalVorticityField(velocities::NamedTuple)

    grid = velocities.u.grid
    computed_dependencies = (velocities.u, velocities.v)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid; computed_dependencies)

    return compute!(Field(ζ_op))
end

function HorizontalFriction(model; ClosureType = AbstractScalarBiharmonicDiffusivity)

    grid          = model.grid
    clock         = model.clock
    closure       = filter(x -> x isa ClosureType, model.closure)
    diffusivities = model.diffusivity_fields
    buoyancy      = model.buoyancy
    velocities    = model.velocities
    free_surface  = model.free_surface
    tracers       = model.tracers
    auxiliary_fields = model.auxiliary_fields

    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    computed_dependencies = (closure, diffusivities, clock, model_fields, buoyancy)

    ∂ⱼ_τ₁ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₁ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₂ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₂ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₃ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₃ⱼ, grid; computed_dependencies)

    τ₁ = compute!(Field(∂ⱼ_τ₁ⱼ_op))
    τ₂ = compute!(Field(∂ⱼ_τ₂ⱼ_op))
    τ₃ = compute!(Field(∂ⱼ_τ₃ⱼ_op))

    return (; τ₁, τ₂, τ₃)
end

function HorizontalFriction(fields::NamedTuple, closure)

    grid          = fields.u.grid
    clock         = Clock{eltype(grid)}(0, 0, 1)
    diffusivities = nothing
    buoyancy      = BuoyancyTracer()
    velocities    = (; u = fields.u, v = fields.v, w = ZFaceField(grid))
    free_surface  = nothing
    tracers       = (; b = fields.b)
    auxiliary_fields = NamedTuple()

    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    computed_dependencies = (closure, diffusivities, clock, model_fields, buoyancy)

    ∂ⱼ_τ₁ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₁ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₂ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₂ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₃ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₃ⱼ, grid; computed_dependencies)

    τ₁ = compute!(Field(∂ⱼ_τ₁ⱼ_op))
    τ₂ = compute!(Field(∂ⱼ_τ₂ⱼ_op))
    τ₃ = compute!(Field(∂ⱼ_τ₃ⱼ_op))

    return (; τ₁, τ₂, τ₃)
end
