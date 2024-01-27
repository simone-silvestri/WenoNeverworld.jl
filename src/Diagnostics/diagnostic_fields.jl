using Oceananigans.Utils
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: immersed_cell
using Oceananigans.Models.HydrostaticFreeSurfaceModels: hydrostatic_fields
using Oceananigans.Coriolis: fᶠᶠᵃ

import Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

##### 
##### Usefull diagnostics
#####

"""
    VerticalVorticity(f::Dict, i)

Returns the three-dimensional vertical vorticity at time index i.
"""
VerticalVorticity(f::Dict, i; indices = (:, :, :)) = compute!(Field(VerticalVorticityOperation(f, i); indices))


"""
    KineticEnergy(f::Dict, i)

Returns the three-dimensional kinetic energy at time index i.
"""
KineticEnergy(f::Dict, i; indices = (:, :, :)) = compute!(Field(KineticEnergyOperation(f, i); indices))


"""
    Stratification(f::Dict, i)

Returns the three-dimensional stratification at time index i.
"""
Stratification(f::Dict, i; indices = (:, :, :)) = compute!(Field(StratificationOperation(f, i); indices))

"""
    PotentialVorticity(f::Dict, i)

Returns the three-dimensional potential vorticity at time index i.
"""
PotentialVorticity(f::Dict, i; indices = (:, :, :)) = compute!(Field(PotentialVorticityOperation(f, i); indices)) 

"""
    DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655)

Returns the three-dimensional density given a buoyancy field b.
"""
DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655, indices = (:, :, :)) = compute!(Field(DensityOperation(b; ρₒ, g); indices))

"""
    DeformationRadius(f::Dict, i)

Returns the two-dimensional deformation vorticity at time index i.
"""
function DeformationRadius(f::Dict, i)
    grid = f[:b].grid
    arch = architecture(grid)

    Ld = Field{Center, Center, Nothing}(grid) 

    launch!(arch, grid, :xy, _deformation_radius!, Ld, f[:b][i], grid, Val(grid.Nz))

    return Ld
end

"""
    VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) 

Returns a three-dimensional field containing the cell volumes at location `loc` with indices `indices`.
"""
VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  
"""
    AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3))

Returns a two-dimensional field containing the cell horizontal areas at location `loc` with indices `indices`.
"""
AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

"""
    HeightField(grid, loc = (Center, Center, Center))  

Returns a three-dimensional field containing the cell vertical spacing at location `loc`.
"""
function HeightField(grid, loc = (Center, Center, Center))  

    zf = Field(loc, grid)
    Lz = grid.Lz

    for k in 1:size(zf, 3)
        interior(zf, :, :, k) .= Lz + znode(k, grid, loc[3]())
    end

    return zf
end

#####
##### KernelFunctionOperations
#####

VerticalVorticityOperation(fields::Dict, i)   =   VerticalVorticityOperation((; u = fields[:u][i], v = fields[:v][i]))
PotentialVorticityOperation(fields::Dict, i)  =  PotentialVorticityOperation((; u = fields[:u][i], v = fields[:v][i], b = fields[:b][i]))
KineticEnergyOperation(fields::Dict, i)       =       KineticEnergyOperation((; u = fields[:u][i], v = fields[:v][i]))
StratificationOperation(fields::Dict, i)      =      StratificationOperation(fields[:b][i])

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

@inline _density_operation(i, j, k, grid, b, ρ₀, g) = ρ₀ * (1 - b[i, j, k] / g)

DensityOperation(b; ρ₀ = 1000.0, g = 9.80655) = 
    KernelFunctionOperation{Center, Center, Center}(_density_operation, b.grid, b, ρ₀, g)

function VerticalVorticityOperation(velocities::NamedTuple)

    grid = velocities.u.grid
    computed_dependencies = (velocities.u, velocities.v)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, computed_dependencies...)

    return ζ_op
end

function StratificationOperation(b)
    grid = b.grid

    N2_op = KernelFunctionOperation{Center, Center, Face}(N²ᶜᶜᶠ, grid, b)

    return N2_op
end

@inline ∂z_bᶠᶠᶜ(i, j, k, grid, b) =  ℑxyzᶠᶠᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
@inline pvᶠᶠᶜ(i, j, k, grid, u, v, b) = (ζ₃ᶠᶠᶜ(i, j, k, grid, u, v) + fᶠᶠᵃ(i, j, k, grid, HydrostaticSphericalCoriolis())) * ∂z_bᶠᶠᶜ(i, j, k, grid, b)

function PotentialVorticityOperation(fields::NamedTuple)

    grid = fields.u.grid
    computed_dependencies = (fields.u, fields.v, fields.b)

    ζ_op = KernelFunctionOperation{Face, Face, Center}(pvᶠᶠᶜ, grid, computed_dependencies...)
    ρ = DensityOperation(fields.b)

    return ζ_op / ρ
end

function KineticEnergyOperation(velocities::NamedTuple)
    u = velocities.u
    v = velocities.v

    E_op = @at (Center, Center, Center) 0.5 * (u^2 + v^2)

    return E_op
end

@inline _deformation_radius(i, j, k, grid, b) = sqrt(max(0, ∂zᶜᶜᶠ(i, j, k, grid, b))) / π /
                                                abs(ℑxyᶜᶜᵃ(i, j, k, grid, fᶠᶠᵃ, HydrostaticSphericalCoriolis()))

@kernel function _deformation_radius!(Ld, b, grid, ::Val{Nz}) where Nz
    i, j = @index(Global, NTuple)

    @inbounds Ld[i, j, 1] = 0

    d₁ᶜᶜᶠ = _deformation_radius(i, j, 1, grid, b)
    @unroll for k in 1:Nz
        d₂ᶜᶜᶠ = _deformation_radius(i, j, k+1, grid, b)
        @inbounds Ld[i, j, k] += ifelse(immersed_cell(i, j, k, grid), 0, 0.5 * (d₁ᶜᶜᶠ + d₂ᶜᶜᶠ) * Δzᶜᶜᶜ(i, j, k, grid))
    end
end
