#####
##### Boundary conditions
#####

using Oceananigans.Utils
using Oceananigans.Grids: node, halo_size
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Operators: Δx, Δy, Az 
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation, HorizontalDivergenceScalarBiharmonicDiffusivity
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.Advection: VorticityStencil, VelocityStencil
using Oceananigans.MultiRegion: multi_region_object_from_array, reconstruct_global_grid

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

@inline surface_wind_stress(i, j, grid, clock, fields, τ) = τ[j]

@inline function grid_specific_array(wind_stress, grid; scaling = 1000.0)

    Ny   = size(grid, 2)
    arch = architecture(grid)
    
    φ_grid = grid.φᵃᶜᵃ[1:Ny]

    τw = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)
        τw[j] = wind_stress(φ, 0.0) ./ scaling
    end

    return arch_array(arch, -τw)
end

@inline function buoyancy_top_relaxation(i, j, grid, clock, fields, p) 

    b = fields.b[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    return @inbounds p.λ * (b - p.initial_buoyancy(x, y, z))
end

@inline function temperature_top_relaxation(i, j, grid, clock, fields, p) 

    T  = fields.T[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)
    Trestoring = p.Sʳ(x, y, z)

    return @inbounds p.λ * (T - Trestoring)
end

@inline function salinity_top_relaxation(i, j, grid, clock, fields, p) 

    S  = fields.S[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    Srestoring = p.Sʳ(x, y, z)
    Sflux      = p.F(x, y, z)

    return @inbounds p.λ * (S - Srestoring) - Sflux
end

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

@inline Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz) =  (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ
@inline    cosine_νhb(i, j, k, grid, lx, ly, lz, clock, fields, ν) = ν * hack_cosd(ynode(ly, j, grid))^3

default_convective_adjustment = RiBasedVerticalDiffusivity()
default_biharmonic_viscosity  = HorizontalDivergenceScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)
default_vertical_diffusivity  = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=1e-5)
default_slope_limiter         = FluxTapering(1e-2)

@inline initialize_model!(model, ::Val{false}, initial_buoyancy, grid, orig_grid, init_file, ::BuoyancyTracer) = set!(model, b = initial_buoyancy)

@inline function initialize_model!(model, ::Val{true}, initial_buoyancy, grid, orig_grid, init_file, ::BuoyancyTracer)
    Hx, Hy, Hz = halo_size(orig_grid)

    b_init = jldopen(init_file)["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    @info "interpolating fields"
    b_init = interpolate_per_level(b_init, orig_grid, grid, (Center, Center, Center))
    u_init = interpolate_per_level(u_init, orig_grid, grid, (Face, Center, Center))
    v_init = interpolate_per_level(v_init, orig_grid, grid, (Center, Face, Center))
    w_init = interpolate_per_level(w_init, orig_grid, grid, (Center, Center, Face))

    set!(model, b=b_init, u=u_init, v=v_init, w=w_init) 
end

@inline initialize_model!(model, ::Val{false}, initial_profiles, grid, orig_grid, init_file, ::SeawaterBuoyancy) = set!(model, T = initial_profiles[1],  S = initial_profiles[2])

@inline function initialize_model!(model, ::Val{true}, initial_temperature, grid, orig_grid, init_file, ::SeawaterBuoyancy)
    Hx, Hy, Hz = halo_size(orig_grid)

    T_init = jldopen(init_file)["T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    S_init = jldopen(init_file)["S/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    @info "interpolating fields"
    T_init = interpolate_per_level(T_init, orig_grid, grid, (Center, Center, Center))
    S_init = interpolate_per_level(S_init, orig_grid, grid, (Center, Center, Center))
    u_init = interpolate_per_level(u_init, orig_grid, grid, (Face, Center, Center))
    v_init = interpolate_per_level(v_init, orig_grid, grid, (Center, Face, Center))
    w_init = interpolate_per_level(w_init, orig_grid, grid, (Center, Center, Face))

    set!(model, b=b_init, u=u_init, v=v_init, w=w_init) 
end

function weno_neverworld_simulation(; grid, 
                                      orig_grid = grid,
                                      μ_drag = 0.001,  
                                      λ_buoy = 7days,
                                      convective_adjustment = default_convective_adjustment,
                                      biharmonic_viscosity  = default_biharmonic_viscosity,
                                      vertical_diffusivity  = default_vertical_diffusivity,
                                      gm_redi_diffusivities = (1000, 1000),
                                      tapering = default_slope_limiter,
                                      coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme()),
                                      free_surface = ImplicitFreeSurface(),
                                      momentum_advection = WENO(vector_invariant = VelocityStencil()),
				                      tracer_advection   = WENO(grid.underlying_grid), 
                                      interp_init = false,
                                      init_file = nothing,
                                      Δt = 5minutes,
                                      stop_time = 10years,
                                      buoyancy_boundary_conditions = true,
                                      velocity_boundary_conditions = true,
                                      initial_buoyancy = initial_buoyancy_tangent,
				                      wind_stress  = zonal_wind_stress,
                                      tracers = :b
                                      )

    # Initializing boundary conditions

    @info "specifying boundary conditions..."

    @apply_regionally τw = grid_specific_array(wind_stress, grid)

    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τw)

    # Quadratic bottom drag:

    if velocity_boundary_conditions
        drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)
        drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)

        u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u) 
        v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v) 

        u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ_drag)
        v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ_drag)

        u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
        v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)
    else
        u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
        v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    end

    Δz_top = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    v_pump = Δz_top / λ_buoy

    b_top_relaxation_bc = FluxBoundaryCondition(buoyancy_top_relaxation, discrete_form=true, parameters = (; λ = v_pump, initial_buoyancy))

    if buoyancy_boundary_conditions
        b_bcs = FieldBoundaryConditions(top = b_top_relaxation_bc)
    else
        b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    end

    #####
    ##### Closures
    #####

    @info "specifying closures..."

    if gm_redi_diffusivities isa Nothing
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)
    else
        κᴳ, κᴿ = gm_redi_diffusivities
        isopycnal_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = κᴳ, κ_symmetric = κᴿ, slope_limiter = tapering)
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment, isopycnal_closure)
    end

    #####
    ##### Model setup
    #####

    @info "building model..."            

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, closure, tracers, momentum_advection, tracer_advection, 
                                          boundary_conditions = (; u = u_bcs, v = v_bcs, b = b_bcs), 
                                          buoyancy = BuoyancyTracer())

    #####
    ##### Model initialization
    #####

    @info "initializing prognostic variables from $(interp_init ? init_file : "scratch")"
    initialize_model!(model, Val(interp_init), initial_buoyancy, grid, orig_grid, init_file, buoyancy_model)

    simulation = Simulation(model; Δt, stop_time)

    @show start_time = [time_ns()]

    function progress(sim)
        sim.model.clock.iteration == 1

        wall_time = (time_ns() - start_time[1]) * 1e-9

        u, v, w = sim.model.velocities

        @info @sprintf("Time: % 12s, it: %d, max(|u|, |v|, |w|): (%.2e, %.2e , %.2e) ms⁻¹, Δt: %.2e s, wall time: %s", 
            prettytime(sim.model.clock.time),
	    sim.model.clock.iteration, maximum(abs, u), maximum(abs, v), maximum(abs, w), sim.Δt,
            prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

    return simulation
end

function neverworld_simulation_seawater(; grid, 
                                          orig_grid = grid,
                                          μ_drag = 0.001,  
                                          λ_T =  7days,
                                          λ_S = 60days,
                                          convective_adjustment = default_convective_adjustment,
                                          biharmonic_viscosity  = default_biharmonic_viscosity,
                                          vertical_diffusivity  = default_vertical_diffusivity,
                                          gm_redi_diffusivities = (1000, 1000),
                                          tapering = default_slope_limiter,
                                          coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme()),
                                          free_surface = ImplicitFreeSurface(),
                                          momentum_advection = WENO(vector_invariant = VelocityStencil()),
                                          tracer_advection   = WENO(grid.underlying_grid), 
                                          interp_init = false,
                                          init_file = nothing,
                                          Δt = 5minutes,
                                          stop_time = 10years,
                                          initial_temperature = initial_temperature_parabola,
                                          initial_salinity = initial_salinity,
                                          salinity_flux = salinity_flux,
                                          equation_of_state = LinearEquationOfState(),
                                          wind_stress  = zonal_wind_stress,
                                          tracers = (:T, :S)
                                          )

    # Initializing boundary conditions

    @info "specifying boundary conditions..."

    @apply_regionally τw = grid_specific_array(wind_stress, grid)

    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τw)

    # Quadratic bottom drag:

    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u) 
    v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v) 

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ_drag)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ_drag)

    u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)

    Δz_top   = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    v_pump_T = Δz_top / λ_T
    v_pump_S = Δz_top / λ_S

    T_top_relaxation_bc = FluxBoundaryCondition(temperature_top_relaxation, discrete_form=true, parameters = (; λ = v_pump_T, Tʳ = initial_temperature))
    S_top_relaxation_bc = FluxBoundaryCondition(salinity_top_relaxation,    discrete_form=true, parameters = (; λ = v_pump_S, Sʳ = initial_salinity, F = salinity_flux))

    T_bcs = FieldBoundaryConditions(top = T_top_relaxation_bc)
    S_bcs = FieldBoundaryConditions(top = S_top_relaxation_bc)
    
    #####
    ##### Closures
    #####

    @info "specifying closures..."

    if gm_redi_diffusivities isa Nothing
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)
    else
        κᴳ, κᴿ = gm_redi_diffusivities
        isopycnal_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = κᴳ, κ_symmetric = κᴿ, slope_limiter = tapering)
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment, isopycnal_closure)
    end

    #####
    ##### Model setup
    #####

    @info "building model..."            

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, closure, tracers, momentum_advection, tracer_advection, 
                                          boundary_conditions = (; u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs), 
                                          buoyancy = SeawaterBuoyancy(; equation_of_state))

    #####
    ##### Model initialization
    #####

    @info "initializing prognostic variables from $(interp_init ? init_file : "scratch")"
    initialize_model!(model, Val(interp_init), (initial_temperature, initial_salinity), grid, orig_grid, init_file, SeawaterBuoyancy())

    simulation = Simulation(model; Δt, stop_time)

    @show start_time = [time_ns()]

    function progress(sim)
        sim.model.clock.iteration == 1

        wall_time = (time_ns() - start_time[1]) * 1e-9

        u, v, w = sim.model.velocities

        @info @sprintf("Time: % 12s, it: %d, max(|u|, |v|, |w|): (%.2e, %.2e , %.2e) ms⁻¹, Δt: %.2e s, wall time: %s", 
        prettytime(sim.model.clock.time),
        sim.model.clock.iteration, maximum(abs, u), maximum(abs, v), maximum(abs, w), sim.Δt,
        prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

    return simulation
end

function run_simulation!(simulation; interp_init = false, init_file = nothing) 
    
    init = interp_init ? true : (init_file isa Nothing ? true : false)

    Δt    = simulation.Δt
    model = simulation.model 
        
    if init
        @info "running simulation from zero-velocity initial conditions"
        run!(simulation)
    else
        @info "running simulation from $init_file"
        update_simulation_clock!(simulation, init_file)
        run!(simulation, pickup=init_file)
    end
    
    @info """
        Simulation took $(prettytime(simulation.run_wall_time))
        Free surface: $(typeof(model.free_surface).name.wrapper)
        Time step: $(prettytime(Δt))
    """
end
