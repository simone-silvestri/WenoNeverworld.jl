using WenoNeverworld
using WenoNeverworld.Auxiliaries
using WenoNeverworld.NeverworldBoundaries
using WenoNeverworld.Parameterizations
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnode
using Test

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_faces(Nz, Depth; h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - Depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
end

@testset "Neverworld Grid" begin
    @info "Testing the Neverworld grid..."
    grid = NeverworldGrid(2)

    @test grid isa Oceananigans.ImmersedBoundaryGrid

    @test grid.Δλᶠᵃᵃ == 2
    @test grid.Δφᵃᶠᵃ == 2
end

@testset "Neverworld Simulation" begin
    @info "Testing the Neverworld simulation..."
    z_faces    = exponential_faces(2, 4000)
    grid       = NeverworldGrid(12; z_faces)
    simulation = weno_neverworld_simulation(grid; stop_iteration = 1)
    run_simulation!(simulation)
end

@inline function c_boundary_condition(i, j, grid, clock, fields) 
    ϕ = φnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    return 1 / 7days * (fields.c[i, j, grid.Nz] - cos(2π * ϕ / grid.Ly))
end

@testset "Tracer Boundary Conditions" begin
    @info "Testing custom tracer boundary conditions..."
    z_faces = exponential_faces(2, 4000)
    grid    = NeverworldGrid(12; z_faces)
    tracers = (:b, :c)

    tracer_boundary_conditions = (; c = c_boundary_condition)

    simulation = weno_neverworld_simulation(grid; stop_iteration = 1, tracers, tracer_boundary_conditions)

    @test simulation.model.tracers.b.boundary_conditions.top.condition.func isa BuoyancyRelaxationBoundaryCondition
    @test simulation.model.tracers.c.boundary_conditions.top.condition.func == c_boundary_condition

    run!(simulation)
end

@testset "Interpolation tests" begin
    @info "Testing three dimensional interpolation and restart from a different grid..."
    # Coarse simulation
    coarse_z_faces = exponential_faces(2, 4000)
    coarse_grid = NeverworldGrid(12; z_faces = coarse_z_faces, H = 1)

    coarse_simulation = weno_neverworld_simulation(coarse_grid; stop_iteration = 1, 
                                                   momentum_advection = nothing, 
                                                   tracer_advection = nothing)
    checkpoint_outputs!(coarse_simulation, "test_fields")

    run!(coarse_simulation)

    b_coarse = coarse_simulation.model.tracers.b

    # Fine simulation interpolated from the coarse one
    fine_z_faces = exponential_faces(4, 4000)
    fine_grid = NeverworldGrid(8; z_faces = fine_z_faces, H = 1)

    @info "    Testing 3-dimensional interpolation..."
    b_fine = regrid_field(b_coarse, coarse_grid, fine_grid, (Center, Center, Center))

    @info "    Testing interpolated restart capabilities..."
    fine_simulation = weno_neverworld_simulation(fine_grid; 
                                                 previous_grid = coarse_grid,
                                                 stop_iteration = 1,
                                                 momentum_advection = nothing, 
                                                 tracer_advection = nothing,
                                                 init_file = "test_fields_checkpoint_iteration0.jld2")

    run!(fine_simulation)
end

@testset "Parameterizations" begin
    @info "Testing parameterization..."

    grid = NeverworldGrid(12; z_faces = [-4000, -2000, 0])
    horizontal_closures = (QGLeith(), EnergyBackScattering())
    for horizontal_closure in horizontal_closures
        @info "    Testing $(typeof(horizontal_closure).name.wrapper) parameterization..."
        simulation = weno_neverworld_simulation(grid; stop_iteration = 1, horizontal_closure)
        run!(simulation)
    end
end