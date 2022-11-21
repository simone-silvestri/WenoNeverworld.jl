using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization

@inline wind_function(y, mid_wind) = 0.0

vertical_diffusivity  = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = 100.0) 
# vertical_diffusivity  = VerticalScalarDiffusivity(ExplicitTimeDiscretization(),           ν = 100.0)
free_surface          = nothing 
momentum_advection    = nothing
tracer_advection      = nothing
biharmonic_viscosity  = nothing
convective_adjustment = nothing
coriolis              = nothing #FPlane(f = -1e-4)

cpugrid = NeverworldGrid(CPU(), 2)
simcpu  = weno_neverworld_simulation(; grid = cpugrid, free_surface, momentum_advection, tracer_advection, coriolis, 
				     wind_stress = wind_function, vertical_diffusivity, biharmonic_viscosity, convective_adjustment)
gpugrid = NeverworldGrid(GPU(), 2)
simgpu  = weno_neverworld_simulation(; grid = gpugrid, free_surface, momentum_advection, tracer_advection, coriolis,
				     wind_stress = wind_function, vertical_diffusivity, biharmonic_viscosity,convective_adjustment)

modcpu = simcpu.model
modgpu = simgpu.model

uᵢ = rand(size(cpugrid))

set!(modcpu, u = uᵢ, b = 0.0)
set!(modgpu, u = uᵢ, b = 0.0)

time_step!(modcpu, 2000)
time_step!(modgpu, 2000)

using Oceananigans.TimeSteppers: update_state!, calculate_tendencies!, ab2_step!, correct_velocities_and_store_tendecies! 

import Oceananigans.Fields: interior

interior(::Nothing) = []

G1cpu = modcpu.timestepper.G⁻
Gncpu = modcpu.timestepper.Gⁿ
G1gpu = modgpu.timestepper.G⁻
Gngpu = modgpu.timestepper.Gⁿ

varcpu = tuple(modcpu.velocities..., modcpu.tracers..., modcpu.pressure.pHY′, G1cpu..., Gncpu...);
vargpu = tuple(modgpu.velocities..., modgpu.tracers..., modgpu.pressure.pHY′, G1gpu..., Gngpu...);

varkeys = tuple(keys(modcpu.velocities)..., :b, :η, :p, keys(G1cpu)..., keys(Gncpu)...)
