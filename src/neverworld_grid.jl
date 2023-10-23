using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size
using Oceananigans.DistributedComputations

"""
    function z_faces_exp(; Nz = 69, Lz = 4000.0, e_folding = 0.06704463421863584)

generates an array of exponential z faces 

"""
function z_faces_exp(; Nz = 69, Lz = 4000.0, e_folding = 0.06704463421863584)
    z_faces   = zeros(Nz + 1)
    Nconstant = 11

    z_faces[1:Nconstant] .= 0:5:50

    for i in 1:(Nz + 1 - Nconstant)
        z_faces[i + Nconstant] = z_faces[i - 1 + Nconstant] + 5 * exp(e_folding * i)
    end

    z_faces    = - reverse(z_faces)
    z_faces[1] = - Lz

    return z_faces
end

"""
    function NeverworldGrid(arch, degree, FT::DataType = Float64; H = 7, longitude = (-2, 62), latitude = (-70, 0), bathymetry_params = NeverWorldBathymetryParameters(), longitudinal_extent = 60) 

builds a `LatitudeLongitudeGrid` with a specified `bathymetry`

Arguments
=========

- `arch` : architecture of the grid, can be `CPU()` or `GPU()` or `Distributed`
- `resolution` : resolution in degrees.
- `FT` : (optional) floating point precision (default = `Float64`)

Keyword Arguments
=================

- `H` : halo size, `Int`
- `longitudinal_extent` : size of the actual domain in longitudinal direction, `Number`
- `longitude` : longitudinal extremes of the domain, `Tuple`. Note: this keyword must be at least `longitude_extent + resolution * 2H`
                to allow for correct advection stencils 
- `latitude` : latitudinal extremes of the domain
- `bathymetry_params` : parameters for the neverworld bathymetry, see `neverworld_bathymetry.jl`
- `z_faces` : array containing the z faces

"""
function NeverworldGrid(resolution, FT::DataType = Float64; 
                        arch = CPU(), H = 7, 
                        longitudinal_extent = 60, 
                        longitude = (-2, 62), 
                        latitude = (-70, 70), 
                        bathymetry_params = NeverWorldBathymetryParameters(),
                        z_faces = z_faces_exp()) 

    Nx = ceil(Int, (longitude[2] - longitude[1]) / resolution)
    Ny = ceil(Int, ( latitude[2] -  latitude[1]) / resolution)
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
        bathy[i, j] = neverworld_bathymetry(λ, φ, bathymetry_params; longitudinal_extent, latitude)
    end

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathy))
end
