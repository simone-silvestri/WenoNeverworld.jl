using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size

const Nz = 69
const Lz = 4000.0

z_faces   = zeros(Nz + 1)
Nconstant = 11

z_faces[1:Nconstant] .= 0:5:50

""" Magic constant """
const t = 0.06704463421863584

for i in 1:(Nz + 1 - Nconstant)
    z_faces[i + Nconstant] = z_faces[i - 1 + Nconstant] + 5 * exp(t * i)
end

z_faces    = - reverse(z_faces)
z_faces[1] = - Lz

function NeverworldGrid(arch, degree, FT::DataType = Float64; H = 5, longitude = (-2, 62), latitude = (-70, 0), bathymetry = bathymetry_without_ridge, longitudinal_extent = 60) 

    Nx = Int((longitude[2] - longitude[1]) / degree)
    Ny = Int((latitude[2]  - latitude[1]) / degree)
    Nz = length(z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny, Nz),
                                            latitude,
                                            longitude,
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = z_faces)

    λ_grid = underlying_grid.λᶜᵃᵃ[1:Nx]
    φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

    bathy = zeros(Nx, Ny)
    for (i, λ) in enumerate(λ_grid), (j, φ) in enumerate(φ_grid)
        bathy[i, j] = bathymetry(λ, φ; longitudinal_extent, latitude)
    end

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathy))
end

