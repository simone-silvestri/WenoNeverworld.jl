using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: VelocityStencil
using UpwindVectorInvariantSchemes
# include("../src/vector_invariant_divergence.jl")

grid = LatitudeLongitudeGrid(size = (20, 10, 1), latitude = (-50, 50), longitude = (-180, 180), z = (-1000, 0), halo = (4, 4, 4))

upwind_scheme      = WENO(VorticityStencil())
momentum_advection = GlobalVectorInvariant(; upwind_scheme)

model = HydrostaticFreeSurfaceModel(; grid, momentum_advection, buoyancy = nothing, coriolis = nothing, closure = nothing)

set!(model, u = (x, y, z) -> cos(y))

u_init = deepcopy(model.velocities.u)

for t in 1:1000
    time_step!(model, 1hour)
end