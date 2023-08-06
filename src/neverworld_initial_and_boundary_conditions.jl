const Lz   = 4000
const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 

""" utility profiles (atan, exponential, and parabolic) """
@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )
@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

struct BuoyancyRelaxationBoundaryCondition{T, S, F} <: Function
    ΔB::T
    λ::S
    func::F
end

BuoyancyRelaxationBoundaryCondition(; ΔB = ΔB, λ = 7days, func = (y, t) -> parabolic_scaling(y)) = BuoyancyRelaxationBoundaryCondition(ΔB, λ, func)

function (b::BuoyancyRelaxationBoundaryCondition)(i, j, grid, clock, fields)
    φ = φnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    b_surf = fields.b[i, j, grid.Nz]
    return 1 / b.λ * (b_surf - b.ΔB * b.func(φ, clock.time))
end

Adapt.adapt_structure(to, b::BuoyancyRelaxationBoundaryCondition) = BuoyancyRelaxationBoundaryCondition(b.ΔB, b.λ, b.func)

struct WindStressBoundaryCondition{F, T, S} <: Function
    φs :: F
    τs :: T
    stress :: S
end

default_φs = (-70, -45, -15, 0, 15, 45, 70)
default_τs = (0.0, 0.2, -0.1, -0.02, -0.1, 0.1, 0.0)
    
WindStressBoundaryCondition(; φs = default_φs, τs = default_τs) =  WindStressBoundaryCondition(φs, τs, nothing)

(ws::WindStressBoundaryCondition)(i, j, grid, clock, fields) = ws.stress[j]

Adapt.adapt_structure(to, ws::WindStressBoundaryCondition) = WindStressBoundaryCondition(nothing, nothing, adapt(to, ws.stress))

"""
    function zonal_wind_stress(y, mid_wind)

returns the zonal wind as per https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf
as a function of latitude `y`
"""
# Fallback!
@inline regularize_boundary_condition(bc,        grid) = bc
@inline regularize_boundary_condition(::Nothing, grid) = zerofunc

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

    return WindStressBoundaryCondition(bc.φs, bc.τs, arch_array(arch, - stress))
end

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

    # Quadratic bottom drag
    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u)
    v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v)

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ_drag)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ_drag)

    u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)

    # Buoyancy boundary conditions
    buoyancy_boundary_condition = regularize_boundary_condition(buoyancy_boundary_condition, grid)
    b_relaxation_bc             = FluxBoundaryCondition(buoyancy_boundary_condition, discrete_form=true)
    b_bcs                       = FieldBoundaryConditions(top = b_relaxation_bc)

    # Additional tracers (outside b)
    tracers = tuple(tracers)
    tracers = filter(tracer -> tracer != :b, tracers)
    tracer_boundary_conditions = validate_tracer_boundary_conditions(tracers, tracer_boundary_conditions)
    tracer_boundary_conditions = materialize_tracer_boundary_conditions(tracers, grid, tracer_boundary_conditions)

    return merge((u = u_bcs, v = v_bcs, b = b_bcs), tracer_boundary_conditions)
end

@inline zerofunc(args...) = 0

function validate_tracer_boundary_conditions(tracers, tracer_boundary_conditions)
    for tracer in tracers
        if !(hasproperty(tracer, tracer_boundary_conditions))
            tracer_boundary_conditions = merge(tracer_boundary_conditions, (; tracer => zerofunc))
        end
    end
end

materialize_tracer_boundary_conditions(tracers::NamedTuple{(), Tuple{}}, args...) = NamedTuple() 

function materialize_tracer_boundary_conditions(tracers, grid, tracer_bcs)
    bcs = NamedTuple()
    for t in tracers
        bc = getproperty(tracer_bcs, t)
        bc = regularize_boundary_condition(bc, grid)
        top_bc = FluxBoundaryCondition(bc, discrete_form=true)
        bcs = merge(bcs, (; t => FieldBoundaryConditions(top = top_bc)))
    end

    return bcs
end




    




