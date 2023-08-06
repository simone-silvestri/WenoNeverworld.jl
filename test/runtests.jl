using WenoNeverworld
using Oceananigans
using Test

@testset "Neverworld Grid" begin
    
end

@testset "Boundary Conditions" begin
    
end

@testset "Neverworld Simulation" begin
    grid       = NeverworldGrid(4)
    simulation = weno_neverworld_simulation(grid; stop_iteration = 1)
    run_simulation!(simulation)
end
