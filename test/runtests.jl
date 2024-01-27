using WenoNeverworld
using WenoNeverworld: exponential_z_faces
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnode
using Test

@testset "Neverworld Grid" begin
    
end

@testset "Neverworld Simulation" begin
    grid       = NeverworldGrid(4)
    simulation = weno_neverworld_simulation(grid; stop_iteration = 1)
    run_simulation!(simulation)
end

@testset "Tracer Boundary Conditions" begin
    grid       = NeverworldGrid(4)
    tracers    = (:b, :c)

    @inline function c_boundary_condition(i, j, grid, clock, fields) 
        ϕ = φnode(i, j, grid.Nz, grid, Center(), Center(), Center())
        return 1 / 7days * (fields.c[i, j, grid.Nz] - cos(2π * ϕ / grid.Ly))
    end

    tracer_boundary_conditions = (; c = c_boundary_condition)

    simulation = weno_neverworld_simulation(grid; stop_iteration = 1, tracers, tracer_boundary_conditions)
    run_simulation!(simulation)
end

@testset "Interpolation tests" begin
    
    # Coarse simulation
    coarse_z_faces = exponential_z_faces(Nz = 10)
    coarse_grid = NeverworldGrid(4; z = coarse_z_faces)

    coarse_simulation = weno_neverworld_simulation(coarse_grid; stop_iteration = 1)
    checkpoint_outputs!(coarse_simulation, "test_fields")

    run_simulation!(coarse_simulation)

    b_coarse = coarse_simulation.model.tracers.b

    # Fine simulation interpolated from the coarse one
    fine_z_faces = exponential_z_faces(Nz = 20)
    fine_grid = Neverworld(2; z = fine_z_faces)

    fine_simulation = weno_neverworld_simulation(fine_grid; 
                                                 previous_grid = coarse_grid,
                                                 init_file = "test_fields_checkpoint_iteration0.jld2")

    run_simulation!(fine_simulation)

    b_fine = WenoNeverworld.regridded_field(b_coarse, fine_grid, (Center, Center, Center))
end