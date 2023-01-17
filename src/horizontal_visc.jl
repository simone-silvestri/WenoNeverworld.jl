using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: HorizontalFormulation
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Operators
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, ℑxyᶜᶜᵃ, ζ₃ᶠᶠᶜ, div_xyᶜᶜᶜ
using Oceananigans.Operators: Δx, Δy
using Oceananigans.Operators: ℑxyz
using CUDA: @allowscalar

@inline Dₛ(i, j, k, grid, u, v) = ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
@inline Dₜ(i, j, k, grid, u, v) = ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
@inline Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz) =  (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))

@inline function νhb_smagorinski_final(i, j, k, grid, lx, ly, lz, clock, fields, p)

   location = (lx, ly, lz)
   from_Dₛ = (Center(), Center(), Center()) 
   from_Dₜ = (Face(),   Face(),   Center()) 
	
   δ₁ = ℑxz(i, j, k, grid, from_Dₛ, location, Dₛ, fields.u, fields.v)    
   δ₂ = ℑxz(i, j, k, grid, from_Dₜ, location, Dₛ, fields.u, fields.v)    

   dynamic_visc = max(p.C * sqrt(δ₁^2 + δ₂^2), 1/p.λ)

   return p.Area(i, j, k, grid, lx, ly, lz)^2 * dynamic_visc
end

function smagorinski_viscosity(formulation; Cₛₘ = 2.0, λ = 5days, Area = Δ²ᵃᵃᵃ)

    @show C = (Cₛₘ / π)^2 / 8

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
 
    A = p.Area(i, j, k, grid, lx, ly, lz)
    visc₁ = dynamic_visc * A^2.5
    visc₂ = A^2 / p.λ

    return max(visc₁, visc₂) 
end

function leith_viscosity(formulation; C_vort = 2.0, C_div = 2.0, λ = 10days, Area = Δ²ᵃᵃᵃ)

    @show C₁ = (C_vort / π)^6 
    @show C₂ = (C_div  / π)^6 

    visc = ScalarBiharmonicDiffusivity(formulation; 
                                       ν=νhb_leith_final, discrete_form=true,  
                                       parameters = (; C₁, C₂, λ, Area))

    @show typeof(visc.ν)

    return visc
end

@inline function νhb_leith_laplacian_final(i, j, k, grid, lx, ly, lz, clock, fields, p)
    
    location = (lx, ly, lz)
    from_∂xζ = (Center(), Face(), Center()) 
    from_∂yζ = (Face(), Center(), Center()) 
    from_∂xδ = (Face(), Center(), Center()) 
    from_∂yδ = (Center(), Face(), Center()) 
	
    ∂xζ = ℑxyz(i, j, k, grid, from_∂xζ, location, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxyz(i, j, k, grid, from_∂yζ, location, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂xδ = ℑxyz(i, j, k, grid, from_∂xδ, location, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑxyz(i, j, k, grid, from_∂yδ, location, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
   
    dynamic_visc = sqrt( p.C₁ * (∂xζ^2 + ∂yζ^2) + p.C₂ * (∂xδ^2 + ∂yδ^2) ) 
 
    A = p.Area(i, j, k, grid, lx, ly, lz)

    return dynamic_visc * A^(3/2)
end

function leith_laplacian_viscosity(formulation = HorizontalFormulation(); C_vort = 1.0, C_div = 1.0, Area = Δ²ᵃᵃᵃ)

    @show C₁ = (C_vort / π)^6 
    @show C₂ = (C_div  / π)^6 

    visc = ScalarDiffusivity(formulation; 
                             ν=νhb_leith_laplacian_final, discrete_form=true,  
                             parameters = (; C₁, C₂, Area))

    @show typeof(visc.ν)

    return visc
end

@inline hack_cosd(φ) = cos(π * φ / 180)
@inline hack_sind(φ) = sin(π * φ / 180)

@inline geometric_νhb(i, j, k, grid, lx, ly, lz, clock, fields, λ) = Δ²ᵃᵃᵃ(i, j, k, grid, lx, ly, lz)^2 / λ
@inline    cosine_νhb(i, j, k, grid, lx, ly, lz, clock, fields, ν) = ν * hack_cosd(ynode(ly, j, grid))^3
