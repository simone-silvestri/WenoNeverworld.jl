using Oceananigans.Utils
using Oceananigans.Grids: node, halo_size, interior_parent_indices
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Operators: Δx, Δy, Az 
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving

using WenoNeverworld.Auxiliaires

#####
##### Default parameterizations for the Neverworld simulation
#####

default_convective_adjustment = RiBasedVerticalDiffusivity()
default_vertical_diffusivity  = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=3e-5)

default_momentum_advection(grid) = VectorInvariant(vorticity_scheme = WENO(order = 9), 
                                                    vertical_scheme = WENO(grid))

"""
    function initialize_model!(model, Val(interpolate), initial_buoyancy, grid, previous_grid, init_file, buoyancymodel)

initializes the model according to interpolate or not on a finer/coarser grid `Val(interpolate)`
"""
@inline initialize_model!(model, ::Val{false}, initial_buoyancy, grid, previous_grid, init_file) = set!(model, b = initial_buoyancy)

@inline function initialize_model!(model, ::Val{true}, initial_buoyancy, grid, previous_grid, init_file)
    Hx, Hy, Hz = halo_size(previous_grid)

    b_init = jldopen(init_file)["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    @info "interpolating fields"
    b_init = regridded_field(b_init, previous_grid, grid, (Center, Center, Center))
    u_init = regridded_field(u_init, previous_grid, grid, (Face, Center, Center))
    v_init = regridded_field(v_init, previous_grid, grid, (Center, Face, Center))
    w_init = regridded_field(w_init, previous_grid, grid, (Center, Center, Face))

    set!(model, b=b_init, u=u_init, v=v_init, w=w_init) 
end

"""
    function weno_neverworld_simulation(grid; 
                                        previous_grid = grid,
                                        μ_drag = 0.001,  
                                        convective_adjustment = default_convective_adjustment,
                                        vertical_diffusivity  = default_vertical_diffusivity,
                                        horizontal_closure    = nothing,
                                        coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()),
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
                                        tracer_boundary_condition = NamedTuple(),
                                        tracers = :b
                                        )

returns a simulation object for the Neverworld simulation.

Arguments:
==========

- `grid`: the grid on which the simulation is to be run

Keyword arguments:
===================
 
- `previous_grid`: the grid on which `init_file` has been generated, if we restart from `init_file`
- `μ_drag`: the drag coefficient for the quadratic bottom drag, default: `0.001`
- `convective_adjustment`: the convective adjustment scheme, default: `RiBasedVerticalDiffusivity()`
- `vertical_diffusivity`: the vertical diffusivity scheme, default: `VerticalScalarDiffusivity(ν=1e-4, κ=3e-5)`
- `horizontal_closure`: the horizontal closure scheme, default: `nothing`
- `coriolis`: the coriolis scheme, default: `HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving())`
- `free_surface`: the free surface scheme, default: SplitExplicitFreeSurface(; grid, cfl = 0.75)
- `momentum_advection`: the momentum advection scheme, default: `VectorInvariant(vorticity_scheme = WENO(order = 9), vertical_scheme = WENO(grid))`
- `tracer_advection`: the tracer advection scheme, default: `WENO(grid)`
- `interp_init`: whether to interpolate the initial conditions from `init_file` to `grid`, default: false
- `init_file`: the file from which to read the initial conditions, default: `nothing`
- `Δt`: the time step, default: `5minutes`
- `stop_time`: the time at which to stop the simulation, default: `10years`
- `stop_iteration`: the iteration at which to stop the simulation, default: Inf
- `initial_buoyancy`: the initial buoyancy field in case of `init_file = nothing`, function of `(x, y, z)` default: `initial_buoyancy_parabola`
- `wind_stress`: the wind stress boundary condition, default: `WindStressBoundaryCondition()` (see `src/neverworld_initial_and_boundary_conditions.jl`)
- `buoyancy_relaxation`: the buoyancy relaxation boundary condition, default: `BuoyancyRelaxationBoundaryCondition()` (see `src/neverworld_initial_and_boundary_conditions.jl`)
- `tracer_boundary_condition`: boundary conditions for tracers outside `:b`, default: nothing
- `tracers`: the tracers to be advected, default: `:b`   
"""
function weno_neverworld_simulation(grid; 
                                    previous_grid = grid,
                                    μ_drag = 0.001,  
                                    convective_adjustment = default_convective_adjustment,
                                    vertical_diffusivity  = default_vertical_diffusivity,
                                    horizontal_closure    = nothing,
                                    coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConserving()),
                                    free_surface = SplitExplicitFreeSurface(; grid, cfl = 0.75),
                                    momentum_advection = default_momentum_advection(grid.underlying_grid),
				                    tracer_advection   = WENO(grid.underlying_grid), 
                                    interp_init = false,
                                    init_file = nothing,
                                    Δt = 5minutes,
                                    stop_time = 10years,
                                    stop_iteration = Inf,
                                    initial_buoyancy = initial_buoyancy_parabola,
				                    wind_stress                = WindStressBoundaryCondition(),
                                    buoyancy_relaxation        = BuoyancyRelaxationBoundaryCondition(),
                                    tracer_boundary_conditions = NamedTuple(),
                                    tracers = :b
                                    )

    # Initializing boundary conditions    
    @info "specifying boundary conditions..."
    boundary_conditions = neverworld_boundary_conditions(grid, μ_drag, wind_stress, buoyancy_relaxation, tracers, tracer_boundary_conditions)

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
    initialize_model!(model, Val(interp_init), initial_buoyancy, grid, previous_grid, init_file)

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

