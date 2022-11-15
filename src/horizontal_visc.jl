using Oceananigans.TurbulenceClosures
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.TurbulenceClosures: ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ, ∂ⱼ_τ₃ⱼ, AbstractScalarBiharmonicDiffusivity
using CUDA: @allowscalar

@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz) =  (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

@inline function HorizontalFriction(model; ClosureType = AbstractScalarBiharmonicDiffusivity)

    grid          = model.grid
    clock         = model.clock
    closure       = filter(x -> x isa ClosureType, model.closure)
    diffusivities = model.diffusivity_fields
    buoyancy      = model.buoyancy
    velocities    = model.velocities
    free_surface  = model.free_surface
    tracers       = model.tracers
    auxiliary_fields = model.auxiliary_fields

    model_fields = merge(hydrostatic_fields(velocities, free_surface, tracers), auxiliary_fields)
    computed_dependencies = (closure, diffusivities, clock, model_fields, buoyancy)

    ∂ⱼ_τ₁ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₁ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₂ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₂ⱼ, grid; computed_dependencies)
    ∂ⱼ_τ₃ⱼ_op = KernelFunctionOperation{Face, Center, Center}(∂ⱼ_τ₃ⱼ, grid; computed_dependencies)

    τ₁ = compute!(Field(∂ⱼ_τ₁ⱼ_op))
    τ₂ = compute!(Field(∂ⱼ_τ₂ⱼ_op))
    τ₃ = compute!(Field(∂ⱼ_τ₃ⱼ_op))

    return (; τ₁, τ₂, τ₃)
end

@inline function νhb_smagorinski_final(i, j, k, grid, lx, ly, lz, clock, fields, p)

   location = (lx, ly, lz)
   from_Dₛ = (Center(), Center(), Center()) 
   from_Dₜ = (Face(),   Face(),   Center()) 
	
   δ₁ = ℑxz(i, j, k, grid, from_Dₛ, location, Dₛ, fields.u, fields.v)    
   δ₂ = ℑxz(i, j, k, grid, from_Dₜ, location, Dₛ, fields.u, fields.v)    

   dynamic_visc = max(p.C * sqrt(δ₁^2 + δ₂^2), 1/p.λ)

   return p.Area(i, j, k, grid, lx, ly, lz)^2 * dynamic_visc
end

function smagorinski_viscosity(formulation, grid; Cₛₘ = 2.0, λ = 5days, Area = Δ²ᵃᵃᵃ)

    dx_min = min_Δx(grid.underlying_grid)
    dy_min = min_Δy(grid.underlying_grid)
    dx_max = @allowscalar grid.Δxᶠᶜᵃ[Int(grid.Ny / 2)]
    dy_max = @allowscalar grid.Δxᶠᶜᵃ[Int(grid.Ny / 2)]
    timescale_max = 100days
    timescale_min = 0.2days

    @show C     = (Cₛₘ / π)^2 / 8
    @show min_ν = (1 / (1 / dx_min^2 + 1 / dy_min^2))^2 / timescale_max
    @show max_ν = (1 / (1 / dx_max^2 + 1 / dy_max^2))^2 / timescale_min

    return ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_smagorinski_final, discrete_form=true,  
				                       parameters = (; C, λ, Area))
end


@inline function νhb_leith_final(i, j, k, grid, lx, ly, lz, clock, fields, p)
    
    location = (lx, ly, lz)
    from_∂xζ = (Center(), Face(), Center()) 
    from_∂yζ = (Face(), Center(), Center()) 
    from_∂xδ = (Face(), Center(), Center()) 
    from_∂yδ = (Center(), Face(), Center()) 
	
    ∂xζ = ℑxyz(i, j, k, grid, from_∂xζ, location, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxyz(i, j, k, grid, from_∂yζ, location, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂xδ = ℑxyz(i, j, k, grid, from_∂xδ, location, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑxyz(i, j, k, grid, from_∂yδ, location, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
   
    dynamic_visc = sqrt( p.C₁ * (∂xζ^2 + ∂yζ^2) + p.C₂ * (∂xδ^2 + ∂yδ^2) ) / 8
 
    A = Area(i, j, k, grid, lx, ly, lz)
    visc₁ = dynamic_visc * A^2.5
    visc₂ = A^2 / p.λ

    return max(visc₁, visc₂) 
end

function leith_viscosity(formulation, grid; C_vort = 2.0, C_div = 2.0, λ = 5days, Area = Δ²ᵃᵃᵃ)

    @show C₁ = (C_vort / π)^6 
    @show C₂ = (C_div  / π)^6 

    visc = ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_leith_final, discrete_form=true,  
                                       parameters = (; C₁, C₂, λ, Area))

    @show typeof(visc.ν)

    return visc
end

@inline νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) =
		(1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))^2 / λ

geometric_viscosity(formulation, timescale) = ScalarBiharmonicDiffusivity(formulation, ν=νhb, discrete_form=true, parameters = timescale) 
