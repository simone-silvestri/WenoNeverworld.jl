"""
    struct QGLeith{FT, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}

QGLeith is a struct representing the Leith scalar diffusivity parameterization for quasi-geostrophic models.

## Fields
- `C`: The coefficient for the diffusivity parameterization.
- `min_N²`: The minimum value for the squared buoyancy frequency.
- `isopycnal_tensor`: The isopycnal tensor model used for the diffusivity calculation.
- `slope_limiter`: The slope limiter used for the diffusivity calculation.

## Constructors
- `QGLeith(FT::DataType = Float64; C=FT(1.0), min_N² = FT(1e-20), isopycnal_model=SmallSlopeIsopycnalTensor(), slope_limiter=FluxTapering(1e-2))`: Construct a QGLeith object with optional parameters.

"""
struct QGLeith{FT, M, S} <: AbstractScalarDiffusivity{ExplicitTimeDiscretization, HorizontalFormulation, 2}
    C :: FT
    min_N² :: FT
    Vscale :: FT
    isopycnal_tensor :: M
    slope_limiter :: S
end

QGLeith(FT::DataType=Float64; C=FT(2), min_N²=FT(1e-20), Vscale=FT(1),
        isopycnal_model=SmallSlopeIsopycnalTensor(), 
        slope_limiter=FluxTapering(1e-2)) =
    QGLeith(C, min_N², Vscale, isopycnal_model, slope_limiter) 

DiffusivityFields(grid, tracer_names, bcs, ::QGLeith) = 
                (; νₑ = CenterField(grid),
                   qʸ  = ZFaceField(grid),
                   qˣ  = ZFaceField(grid),
                   Ld = Field{Center, Center, Nothing}(grid))

@inline function ∇h_ζ(i, j, k, grid, coriolis, fields)

    ∂xζ = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)

    ∂xf = ℑyᵃᶜᵃ(i, j, k, grid, ∂xᶜᶠᶜ, fᶠᶠᵃ, coriolis)
    ∂yf = ℑxᶜᵃᵃ(i, j, k, grid, ∂yᶠᶜᶜ, fᶠᶠᵃ, coriolis)
    
    return ∂xζ + ∂xf, ∂yζ + ∂yf
end

@inline function abs²_∇h_δ(i, j, k, grid, fields)

    ∂xδ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    
    return (∂xδ^2 + ∂yδ^2)
end

@kernel function calculate_qgleith_viscosity!(ν, Ld, qˣ, qʸ, grid, closure, velocities, coriolis)
    i, j, k = @index(Global, NTuple)

    ∂ζx, ∂ζy =  ∇h_ζ(i, j, k, grid, coriolis, velocities)
    ∂qx = ∂zᶜᶜᶜ(i, j, k, grid, qˣ)
    ∂qy = ∂zᶜᶜᶜ(i, j, k, grid, qʸ) 

    ∂δ² =  abs²_∇h_δ(i, j, k, grid, velocities)

    fᶜᶜᶜ = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis)

    ∂ζ² = ∂ζx^2 + ∂ζy^2
    ∂q² = (∂qx + ∂ζx)^2 + (∂qy + ∂ζy)^2

    A  = Δ²ᶜᶜᶜ(i, j, k, grid)
    Δs = A^0.5

    Bu  = Ld[i, j, 1]^2 / A
    Ro  = closure.Vscale / (fᶜᶜᶜ * Δs)
    
    ∂Q² = min(∂q², ∂ζ² * (1 + 1 / Bu)^2)
    ∂Q² = min(∂Q², ∂ζ² * (1 + 1 / Ro^2)^2)

    C = closure.C

    @inbounds ν[i, j, k] = (C * Δs / π)^(3) * sqrt(∂Q² + ∂δ²) 
end

@inline ∂yb_times_f2_div_N2(i, j, k, grid, clo, coriolis, buoyancy, tracers) = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) / 
                                                                               max(clo.min_N², ∂z_b(i, j, k, grid, buoyancy, tracers)) *
                                                                               ℑyzᵃᶜᶠ(i, j, k, grid, ∂y_b, buoyancy, tracers)

@inline ∂xb_times_f2_div_N2(i, j, k, grid, clo, coriolis, buoyancy, tracers) = ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) / 
                                                                               max(clo.min_N², ∂z_b(i, j, k, grid, buoyancy, tracers))  *
                                                                               ℑxzᶜᵃᶠ(i, j, k, grid, ∂x_b, buoyancy, tracers)

@kernel function compute_stretching!(qˣ, qʸ, grid, closure, tracers, buoyancy, coriolis)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        qˣ[i, j, k] = ∂xb_times_f2_div_N2(i, j, k, grid, closure, coriolis, buoyancy, tracers)
        qʸ[i, j, k] = ∂yb_times_f2_div_N2(i, j, k, grid, closure, coriolis, buoyancy, tracers)
    end
end

@inline _deformation_radius(i, j, k, grid, C, buoyancy, coriolis) = sqrt(max(0, ∂z_b(i, j, k, grid, buoyancy, C))) / π /
                                                                         abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis))

@kernel function calculate_deformation_radius!(Ld, grid, tracers, buoyancy, coriolis)
    i, j = @index(Global, NTuple)

    @inbounds begin
        Ld[i, j, 1] = 0
        @unroll for k in 1:grid.Nz
            Ld[i, j, 1] += Δzᶜᶜᶠ(i, j, k, grid) * _deformation_radius(i, j, k, grid, tracers, buoyancy, coriolis)
        end
    end
end

function compute_diffusivities!(diffusivity_fields, closure::QGLeith, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    coriolis = model.coriolis

    launch!(arch, grid, :xy, 
            calculate_deformation_radius!, diffusivity_fields.Ld, grid, tracers, buoyancy, coriolis)

    launch!(arch, grid, parameters,
            compute_stretching!, diffusivity_fields.qˣ, diffusivity_fields.qʸ, grid, closure, tracers, buoyancy, coriolis)

    launch!(arch, grid, parameters,
            calculate_qgleith_viscosity!,
            diffusivity_fields.νₑ, diffusivity_fields.Ld, 
            diffusivity_fields.qˣ, diffusivity_fields.qʸ, 
            grid, closure, velocities, coriolis)

    return nothing
end

@inline viscosity(::QGLeith, K) = K.νₑ
@inline diffusivity(::QGLeith, K, ::Val{id}) where id = K.νₑ   

#####
##### Abstract Smagorinsky functionality
#####

@inline diffusive_flux_x(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, args...) where tracer_index = zero(grid)
@inline diffusive_flux_y(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, args...) where tracer_index = zero(grid)
@inline diffusive_flux_z(i, j, k, grid, closure::QGLeith, diffusivities, ::Val{tracer_index}, args...) where tracer_index = zero(grid)

#=
@inline function diffusive_flux_x(i, j, k, grid, closure::QGLeith, diffusivities, 
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ    = diffusivities.νₑ
    νₑⁱʲᵏ = ℑxᶠᵃᵃ(i, j, k, grid, νₑ)
    ∂x_c  = ∂xᶠᶜᶜ(i, j, k, grid, c)
    ∂z_c  = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, c)

    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)
    R₁₃ = isopycnal_rotation_tensor_xz_fcc(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)

    return - νₑⁱʲᵏ * ϵ * (∂x_c + R₁₃ * ∂z_c)  
end

@inline function diffusive_flux_y(i, j, k, grid, closure::QGLeith, diffusivities,
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ    = diffusivities.νₑ
    νₑⁱʲᵏ = ℑyᵃᶠᵃ(i, j, k, grid, νₑ)
    ∂y_c  = ∂yᶜᶠᶜ(i, j, k, grid, c)
    ∂z_c  = ℑyzᵃᶠᶜ(i, j, k, grid, ∂zᶜᶜᶠ, c)

    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)
    R₂₃ = isopycnal_rotation_tensor_yz_cfc(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)

    return - νₑⁱʲᵏ * ϵ * (∂y_c + R₂₃ * ∂z_c)  
end

@inline function diffusive_flux_z(i, j, k, grid, closure::QGLeith, diffusivities, 
                                  ::Val{tracer_index}, c, clock, fields, buoyancy) where tracer_index

    νₑ = diffusivities.νₑ

    νₑⁱʲᵏ = ℑzᵃᵃᶠ(i, j, k, grid, νₑ)

    ∂x_c = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, c)
    ∂y_c = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᶜᶠᶜ, c)
    ∂z_c = ∂zᶜᶜᶠ(i, j, k, grid, c)

    R₃₁ = isopycnal_rotation_tensor_xz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)
    R₃₂ = isopycnal_rotation_tensor_yz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)
    R₃₃ = isopycnal_rotation_tensor_zz_ccf(i, j, k, grid, buoyancy, fields, closure.isopycnal_tensor)

    ϵ = tapering_factor(i, j, k, grid, closure, fields, buoyancy)

    return - νₑⁱʲᵏ * ϵ * (R₃₁ * ∂x_c +
                          R₃₂ * ∂y_c + 
                          R₃₃ * ∂z_c)
end
=#