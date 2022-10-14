#####
##### Boundary conditions
#####

underlying_grid = grid.underlying_grid

Nx, Ny, Nz = size(grid)

λ_grid = underlying_grid.λᶜᵃᵃ[1:Nx]
φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

τw = zeros(Ny)
for (j, φ) in enumerate(φ_grid)
    τw[j] = wind_stress(φ, 0.0) ./ 1000
end
τw = arch_array(arch, -τw)

@inline surface_wind_stress(i, j, grid, clock, fields, τ) = τ[j]
u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τw)

using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ

# Quadratic bottom drag:
μ = 0.003 # Non dimensional

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)^2)^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)^2)^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ)
drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ)

u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u)
v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v)

u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)

Δz_top = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
λ = Δz_top / 7days

@inline buoyancy_top_relaxation(x, y, z, t, b, p) = @inbounds p.λ * (b - p.initial_bouyancy(x, y, z))

b_top_relaxation_bc = FluxBoundaryCondition(buoyancy_top_relaxation, field_dependencies = :b, parameters = (; λ, initial_buoyancy))

b_bcs = FieldBoundaryConditions(top = b_top_relaxation_bc)

#####
##### Closures
#####

using Oceananigans.Operators: Δx, Δy
using Oceananigans.TurbulenceClosures

@inline νhb(i, j, k, grid, lx, ly, lz, clock, fields) = (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2 ))^2 / 5days

vertical_diffusivity  = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)
convective_adjustment = RiBasedVerticalDiffusivity()
biharmonic_viscosity  = HorizontalDivergenceScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true)

closures = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)

#####
##### Model setup
#####

free_surface = ImplicitFreeSurface(solver_method=:HeptadiagonalIterativeSolver)

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    free_surface = free_surface,
                                    momentum_advection = WENO(vector_invariant = VelocityStencil()),
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    closure = closures,
                                    boundary_conditions = (; u = u_bcs, v = v_bcs, b = b_bcs),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    tracer_advection = WENO(underlying_grid))

#####
##### Model initialization
#####

if !interp_init
    set!(model, b = initial_buoyancy)
else
    set!(model, b=b_init, u=u_init, v=v_init, w=w_init, η=η_init) 
end

simulation = Simulation(model; Δt, stop_time)

start_time = [time_ns()]

u, v, w = model.velocities
b = model.tracers.b
η = model.free_surface.η

output_fields = (; u, v, w, b, η)

ke = Field(u^2 + v^2)
u2 = Field(u^2)
v2 = Field(v^2)
η2 = Field(η^2)
b2 = Field(b^2)
w2 = Field(w^2)
vb = Field(v * b)
ub = Field(u * b)
wb = Field(w * b)

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u
    compute!(ke)

    @info @sprintf("Time: % 12s, iteration: %d, max(|u|, |v|, |w|): (%.2e, %.2e , %.2e) ms⁻¹ max(ke): %.2e m²s⁻², Δt: %.2e s, wall time: %s", 
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(ke), sim.Δt,
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard]   = Callback(wizard, IterationInterval(50))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

ζ  = VerticalVorticityField(model)
ζ2 = Field(ζ^2)

averaged_fields = (; u, v, w, b, η, ke, ζ, ζ2, u2, v2, w2, η2, b2, vb, ub, wb)

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, output_fields,
                                                              schedule = TimeInterval(30days),
                                                              filename = output_prefix * "_snapshot",
                                                              overwrite_existing = true)

simulation.output_writers[:averaged_fields] = JLD2OutputWriter(model, averaged_fields,
                                                               schedule = AveragedTimeInterval(30days, window=30days, stride = 4),
                                                               filename = output_prefix * "_averaged",
                                                               overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(1year),
                                                        prefix = output_prefix * "_checkpoint",
                                                        overwrite_existing = true)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

if init
    run!(simulation)
else
    run!(simulation, pickup=init_file)
end

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""