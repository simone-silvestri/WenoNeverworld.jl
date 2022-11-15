#####
##### Boundary conditions
#####

using CUDA

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

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ)
drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ)

no_slip = FluxBoundaryCondition(0.0)

u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u) #, north = no_slip, south = no_slip, west = no_slip, east = no_slip)
v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v) #, north = no_slip, south = no_slip, west = no_slip, east = no_slip)

u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)

Δz_top = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
λ = Δz_top / 7days

using Oceananigans.Grids: node

@inline function buoyancy_top_relaxation(i, j, grid, clock, fields, p) 
    
    b = fields.b[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    return @inbounds p.λ * (b - p.initial_buoyancy(x, y, z))
end

b_top_relaxation_bc = FluxBoundaryCondition(buoyancy_top_relaxation, discrete_form=true, parameters = (; λ, initial_buoyancy))

b_bcs = FieldBoundaryConditions(top = b_top_relaxation_bc)

#####
##### Closures
#####

using Oceananigans.Operators: Δx, Δy, Az
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation, HorizontalDivergenceScalarBiharmonicDiffusivity

include("horizontal_visc.jl")

@inline νhb_old(i, j, k, grid, lx, ly, lz, clock, fields, λ) =
                (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))^2 / λ

@inline νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Az(i, j, k, grid, lx, ly, lz)^2 / λ

biharmonic_viscosity   = HorizontalDivergenceScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters = λ)

vertical_diffusivity  = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-4, κ=1e-5)
convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), convective_κz = 0.2, convective_νz = 0.5)

closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)

#####
##### Model setup
#####

using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.Advection: VorticityStencil, VelocityStencil

coriolis = HydrostaticSphericalCoriolis()

free_surface = ImplicitFreeSurface(solver_method=:HeptadiagonalIterativeSolver)

momentum_advection = WENO(vector_invariant = VelocityStencil())

model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, closure, momentum_advection,
                                    boundary_conditions = (; u = u_bcs, v = v_bcs, b = b_bcs),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    tracer_advection = WENO(underlying_grid))
        
#####
##### Model initialization
#####

if !interp_init
    set!(model, b=initial_buoyancy)
else
    b_init = jldopen(init_file)["b/data"][H+1:end-H, H+1:end-H, H+1:end-H]
    u_init = jldopen(init_file)["u/data"][H+1:end-H, H+1:end-H, H+1:end-H]
    v_init = jldopen(init_file)["v/data"][H+1:end-H, H+1:end-H, H+1:end-H]
    if !(grid == orig_grid)
         @info "interpolating b field"
         b_init = interpolate_per_level(b_init, old_degree, new_degree, (Center, Center, Center), H)
         @info "interpolating u field"
         u_init = interpolate_per_level(u_init, old_degree, new_degree, (Face, Center, Center), H)
         @info "interpolating v field"
         v_init = interpolate_per_level(v_init, old_degree, new_degree, (Center, Face, Center), H)
    end
    set!(model, b=b_init, u=u_init, v=v_init) 
end

simulation = Simulation(model; Δt, stop_time)

start_time = [time_ns()]

u, v, w = model.velocities
b = model.tracers.b
η = model.free_surface.η

output_fields = (; u, v, w, b)

η2 = η^2
u2 = u^2
v2 = v^2
b2 = b^2
w2 = w^2
vb = v * b
ub = u * b
wb = w * b

using Statistics: mean

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u

    @info @sprintf("Time: % 12s, it: %d, max(|u|, |v|, |w|): (%.2e, %.2e , %.2e) ms⁻¹, Δt: %.2e s, wall time: %s", 
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u), maximum(abs, v), maximum(abs, w), sim.Δt,
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

ζ  = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid; computed_dependencies = (u, v))
ζ2 = ζ^2

averaged_fields = (; u, v, w, b, ζ, ζ2, u2, v2, w2, b2, ub, vb, wb)

simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields,
                                                              schedule = TimeInterval(30days),
                                                              filename = output_prefix * "_snapshots",
                                                              overwrite_existing = true)

simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, (u, v, w, b),
                                                              schedule = TimeInterval(5days),
                                                              filename = output_prefix * "_surface",
                                                              indices = (:, :, grid.Nz),
                                                              overwrite_existing = true)

simulation.output_writers[:averaged_fields] = JLD2OutputWriter(model, averaged_fields,
                                                               schedule = AveragedTimeInterval(30days, window=30days, stride = 10),
                                                               filename = output_prefix * "_averages",
                                                               overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(checkpoint_time),
                                                        prefix = output_prefix * "_checkpoint",
                                                        overwrite_existing = true)
