#####
##### Boundary conditions
#####

using Oceananigans.Utils
using Oceananigans.Grids: node, halo_size, interior_parent_indices
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Operators: Δx, Δy, Az 
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation, HorizontalDivergenceScalarBiharmonicDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using Oceananigans.MultiRegion: multi_region_object_from_array, reconstruct_global_grid

#####
##### Default parameterizations for the Neverworld simulation
#####

default_convective_adjustment = RiBasedVerticalDiffusivity()
default_vertical_diffusivity  = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=1e-5)

default_momentum_advection(grid) = VectorInvariant(vorticity_scheme = WENO(order = 9), 
                                                    vertical_scheme = WENO(grid))

"""
    function initialize_model!(model, Val(interpolate), initial_buoyancy, grid, orig_grid, init_file, buoyancymodel)

initializes the model according to interpolate or not on a finer/coarser grid `Val(interpolate)`
"""
@inline initialize_model!(model, ::Val{false}, initial_buoyancy, grid, orig_grid, init_file) = set!(model, b = initial_buoyancy)

@inline function initialize_model!(model, ::Val{true}, initial_buoyancy, grid, orig_grid, init_file)
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

function weno_neverworld_simulation(grid; 
                                    orig_grid = grid,
                                    μ_drag = 0.001,  
                                    convective_adjustment = default_convective_adjustment,
                                    vertical_diffusivity  = default_vertical_diffusivity,
                                    horizontal_closure    = nothing,
                                    coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme()),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl = 0.75),
                                    momentum_advection = default_momentum_advection(grid.underlying_grid),
				                    tracer_advection   = WENO(grid.underlying_grid), 
                                    interp_init = false,
                                    init_file = nothing,
                                    Δt = 5minutes,
                                    stop_time = 10years,
                                    stop_iteration = Inf,
                                    initial_buoyancy = initial_buoyancy_parabola,
				                    wind_stress               = WindStressBoundaryCondition(),
                                    buoyancy_relaxation       = BuoyancyRelaxationBoundaryCondition(),
                                    tracer_boundary_condition = nothing,
                                    tracers = :b
                                    )

    # Initializing boundary conditions

    @info "specifying boundary conditions..."
    boundary_conditions = neverworld_boundary_conditions(grid, μ_drag, wind_stress, buoyancy_relaxation, tracers, tracer_boundary_condition)

    #####
    ##### Closures
    #####

    @info "specifying closures..."
    closure = (vertical_diffusivity, horizontal_closure, convective_adjustment)

    #####
    ##### Model setup
    #####

    @info "building model..."            
    model = HydrostaticFreeSurfaceModel(; grid, free_surface, 
                                          coriolis,
                                          closure, 
                                          tracers, 
                                          momentum_advection, 
                                          tracer_advection, 
                                          boundary_conditions, 
                                          buoyancy = BuoyancyTracer())

    #####
    ##### Model initialization
    #####

    @info "initializing prognostic variables from $(interp_init ? init_file : "scratch")"
    initialize_model!(model, Val(interp_init), initial_buoyancy, grid, orig_grid, init_file)

    simulation = Simulation(model; Δt, stop_time, stop_iteration)

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

