using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, AbstractScalarBiharmonicDiffusivity

import Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

function VerticalVorticityField(velocities::Dict, i)

    grid = velocities[:u].grid
    
    dependencies = (velocities[:u][i], velocities[:v][i])
    fill_halo_regions!(dependencies)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid; computed_dependencies = dependencies)

    return compute!(Field(ζ_op))
end


function KineticEnergyField(fields::Dict, i)
    fill_halo_regions!(fields[:u][i])
    fill_halo_regions!(fields[:v][i])

    E_op = @at (Center, Center, Center) 0.5 * (fields[:u][i]^2 + fields[:v][i]^2)

    return compute!(Field(E_op))
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
