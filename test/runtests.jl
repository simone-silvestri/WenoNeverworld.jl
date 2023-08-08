using WenoNeverworld
using Oceananigans
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
