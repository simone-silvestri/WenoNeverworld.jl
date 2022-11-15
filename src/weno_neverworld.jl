#####
##### Boundary conditions
#####

using Oceananigans.Grids: node, halo_size
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Operators: Δx, Δy, Az 
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation, HorizontalDivergenceScalarBiharmonicDiffusivity
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.Advection: VorticityStencil, VelocityStencil

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

@inline surface_wind_stress(i, j, grid, clock, fields, τ) = τ[j]

@inline function buoyancy_top_relaxation(i, j, grid, clock, fields, p) 

    b = fields.b[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    return @inbounds p.λ * (b - p.initial_buoyancy(x, y, z))
end

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Az(i, j, k, grid, lx, ly, lz)^2 / λ

default_convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 0.2, convective_νz = 0.5)
default_biharmonic_viscosity  = HorizontalDivergenceScalarBiharmonicDiffusivity(ν=geometric_νhb, discrete_form=true, parameters = 5days)
default_vertical_diffusivity  = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)
default_slope_limiter         = FluxTapering(1e-2)

function weno_neverworld_simulation(; grid, orig_grid,
                                      μ_drag = 0.003,  
                                      λ_buoy = 7days,
                                      convective_adjustment = default_convective_adjustment,
                                      biharmonic_viscosity  = default_biharmonic_viscosity,
                                      vertical_diffusivity  = default_vertical_diffusivity,
                                      gm_redi_diffusivities = nothing,
                                      tapering = default_slope_limiter,
                                      coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme()),
                                      momentum_advection = WENO(vector_invariant = VelocityStencil()),
                                      interp_init = false,
                                      init_file = nothing,
                                      Δt = 5minutes,
                                      stop_time = 10years
                                      )

    arch = architecture(grid)
    underlying_grid = grid.underlying_grid

    Nx, Ny, Nz = size(grid)

    φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

    # Initializing boundary conditions

    @info "specifying boundary conditions..."

    τw = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)
        τw[j] = wind_stress(φ, 0.0) ./ 1000
    end
    τw = arch_array(arch, -τw)

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

    Δz_top = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    v_pump = Δz_top / λ_buoy

    b_top_relaxation_bc = FluxBoundaryCondition(buoyancy_top_relaxation, discrete_form=true, parameters = (; λ = v_pump, initial_buoyancy))

    b_bcs = FieldBoundaryConditions(top = b_top_relaxation_bc)

    #####
    ##### Closures
    #####

    @info "specifying closures..."

    if gm_redi_diffusivities isa Nothing
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)
    else
        κᴳ, κᴿ = gm_redi_diffusivities
        isopycnal_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = κᴳ, κ_symmetric = κᴿ, tapering)
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment, isopycnal_closure)
    end

    #####
    ##### Model setup
    #####

    @info "building model..."
    
    free_surface = ImplicitFreeSurface(solver_method=:HeptadiagonalIterativeSolver)

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, closure, momentum_advection,
                                        boundary_conditions = (; u = u_bcs, v = v_bcs, b = b_bcs),
                                        buoyancy = BuoyancyTracer(),
                                        tracers = :b,
                                        tracer_advection = WENO(underlying_grid))

    #####
    ##### Model initialization
    #####

    @info "initializing prognostic variables from $(interp_init ? init_file : "scratch")"

    if !interp_init
        set!(model, b=initial_buoyancy)
    else
        Hx, Hy, Hz = halo_size(orig_grid)
        b_init = jldopen(init_file)["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        if !(grid == orig_grid)
             @info "interpolating b field"
             b_init = interpolate_per_level(b_init, orig_grid, grid, (Center, Center, Center))
             @info "interpolating u field"
             u_init = interpolate_per_level(u_init, orig_grid, grid, (Face, Center, Center))
             @info "interpolating v field"
             v_init = interpolate_per_level(v_init, orig_grid, grid, (Center, Face, Center))
        end
        set!(model, b=b_init, u=u_init, v=v_init) 
    end

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

function run_simulation!(simulation; init = true, init_file = nothing) 
    
    Δt    = simulation.Δt
    model = simulation.model 
        
    if init
        run!(simulation)
    else
        update_simulation_clock!(simulation, init_file)
        run!(simulation, pickup=init_file)
    end
    
    @info """
        Simulation took $(prettytime(simulation.run_wall_time))
        Free surface: $(typeof(model.free_surface).name.wrapper)
        Time step: $(prettytime(Δt))
    """
end