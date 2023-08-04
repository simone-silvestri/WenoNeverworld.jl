const Lz   = 4000
const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const fact = 5.0

""" utility profiles (atan, exponential, and parabolic) """
@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )
@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline atan_scaling(y)      = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2

@inline initial_buoyancy_tangent(x, y, z)  = exponential_profile(z) * atan_scaling(y)
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

@Base.kwdef struct BuoyancyRelaxationBoundaryCondition{T, S}
    ΔB::T
    λ::S
end

BuoyancyRelaxationBoundaryCondition() = BuoyancyRelaxationBoundaryCondition(7 * 86400, 6.0e-2)

function (b::BuoyancyRelaxationBoundaryCondition)(i, j, grid, clock, fields)
    φ = φnode(i, j, grid.Nz, Center(), Center(), Center())
    b_surf = fields.b[i, j, grid.Nz]
    return 1 / b.λ * (b_surf - b.ΔB * parabolic_scaling(φ))
end

@Base.kwdef struct WindStressBoundaryCondition{F, S, T}
    func :: F
    stress :: S
    params :: P
end
    
WindStressBoundaryCondition() = WindStressBoundaryCondition(zonal_wind_stress, default_wind_stress_parameters)

(ws::WindStressBoundaryCondition)(i, j, grid, clock, fields) = ws.stress[j]

"""
    function zonal_wind_stress(y, mid_wind)

returns the zonal wind as per https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf
as a function of latitude `y`
"""

function default_wind_stress_parameters()
    φs = [-45, -15, 0, 15, 45]
    ntlist = NamedTuple{(:x₁, :x₂, :y₁, :y₂),NTuple{4,Float64}}[]
    push!(ntlist, (; x₁ = -70.0, x₂ = -45.0, y₁ = 0.0, y₂ = 0.2))
    push!(ntlist, (; x₁ = -45.0, x₂ = -15.0, y₁ = 0.2, y₂ = -0.1))
    push!(ntlist, (; x₁ = -15.0, x₂ = 0.0, y₁ = -0.1, y₂ = -0.02))
    push!(ntlist, (; x₁ = 0.0, x₂ = 15.0, y₁ = -0.02, y₂ = -0.1))
    push!(ntlist, (; x₁ = 15.0, x₂ = 45.0, y₁ = -0.1, y₂ = 0.1))
    push!(ntlist, (; x₁ = 45.0, x₂ = 70.0, y₁ = 0.1, y₂ = 0.0))
    return (; φs, ntlist)
end

#TODO: ADD CHECK TO MAKE SURE THIS IS CORRECT
@inline function zonal_wind_stress(φ, p)
    #=
    if φ < -45
        return cubic_profile(φ, x₁ = -70.0, x₂ = -45.0, y₁ = 0.0, y₂ = 0.2) / 1000
    elseif φ < -15
        return cubic_profile(φ, x₁ = -45.0, x₂ = -15.0, y₁ = 0.2, y₂ = -0.1) / 1000
    elseif φ < 0
        return cubic_profile(φ, x₁ = -15.0, x₂ = 0.0, y₁ = -0.1, y₂ = -0.02) / 1000
    elseif φ < 15
        return cubic_profile(φ, x₁ = 0.0, x₂ = 15.0, y₁ = -0.02, y₂ = -0.1) / 1000
    elseif φ < 45
        return cubic_profile(φ, x₁ = 15.0, x₂ = 45.0, y₁ = -0.1, y₂ = 0.1) / 1000
    else
        return cubic_profile(φ, x₁ = 45.0, x₂ = 70.0, y₁ = 0.1, y₂ = 0.0) / 1000
    end
    =#
    (; φs, ntlist) = p
    varphi_index = sum(φ .> φs) + 1
    return cubic_profile(φ; ntlist[varphi_index]...) / 1000
end

@inline function regularize_top_boundary_condition(bc::WindStressBoundaryCondition, grid)

    Ny   = size(grid, 2)
    arch = architecture(grid)
    
    φ_grid = grid.φᵃᶜᵃ[1:Ny]

    bc_array = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)
        bc_array[j] = bc(φ, bc.params) 
    end

    return WindStressBoundaryCondition(bc.func, arch_array(arch, - bc_array), bc.params)
end

# Fallback!
@inline regularize_top_boundary_condition(bc_array, grid) = nothing

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

function neverworld_boundary_conditions(grid, μ_drag, wind_stress, buoyancy_boundary_condition, tracers, tracer_boundary_conditions)
    
    # Velocity boundary conditions
    wind_stress       = regularize_top_boundary_condition(wind_stress, grid)
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
    buoyancy_boundary_condition = regularize_top_boundary_condition(buoyancy_boundary_condition, grid)
    b_relaxation_bc             = FluxBoundaryCondition(buoyancy_boundary_condition, discrete_form=true)
    b_bcs                       = FieldBoundaryConditions(top = b_relaxation_bc)

    # Additional tracers (outside b)
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
        bc = regularize_top_boundary_condition(bc, grid)
        bcs = merge(bcs, (; t => FluxBoundaryCondition(bc, discrete_form=true)))
    end

    return bcs
end




    




