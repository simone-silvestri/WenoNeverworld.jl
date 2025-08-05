module NeverworldBoundaries

export neverworld_boundary_conditions
export BuoyancyRelaxationBoundaryCondition
export WindStressBoundaryCondition
export initial_buoyancy_parabola

using WenoNeverworld
using WenoNeverworld.Auxiliaries
using WenoNeverworld.NeverworldGrids
using WenoNeverworld.Constants
using WenoNeverworld.Auxiliaries: parabolic_scaling, exponential_profile

using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: interpolate
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: λnode, φnode, halo_size, on_architecture
using Oceananigans.Utils: instantiate
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using Adapt

# Fallback!
@inline regularize_boundary_condition(bc,        grid) = bc
@inline regularize_boundary_condition(::Nothing, grid) = zerofunc

include("buoyancy_relaxation_bc.jl")
include("wind_stress_bc.jl")
include("tracer_boundary_conditions.jl")

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

function neverworld_boundary_conditions(grid, μ_drag, wind_stress, buoyancy_boundary_condition, tracers, tracer_boundary_conditions)
    
    # Velocity boundary conditions
    wind_stress       = regularize_boundary_condition(wind_stress, grid)
    u_wind_stress_bc  = FluxBoundaryCondition(wind_stress, discrete_form=true)

    if μ_drag > 0
        # Quadratic bottom drag
        drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)
        drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)

        u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u)
        v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v)

        u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ_drag)
        v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ_drag)
        
        u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
        v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)
    else
        u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc)
        v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(nothing))
    end
    
    # Buoyancy boundary conditions
    buoyancy_boundary_condition = regularize_boundary_condition(buoyancy_boundary_condition, grid)
    b_relaxation_bc             = FluxBoundaryCondition(buoyancy_boundary_condition, discrete_form=true)
    b_bcs                       = FieldBoundaryConditions(top = b_relaxation_bc)

    # Additional tracers (outside b)
    tracers = tracers isa Symbol ? tuple(tracers) : tracers
    tracers = filter(tracer -> tracer != :b, tracers)
    tracer_boundary_conditions = validate_tracer_boundary_conditions(tracers, tracer_boundary_conditions)
    tracer_boundary_conditions = materialize_tracer_boundary_conditions(tracers, grid, tracer_boundary_conditions)

    return merge((u = u_bcs, v = v_bcs, b = b_bcs), tracer_boundary_conditions)
end

end


