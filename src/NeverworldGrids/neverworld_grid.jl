using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size
using Oceananigans.DistributedComputations

"""
    function exponential_z_faces(; Nz = 69, Lz = 4000.0, e_folding = 0.06704463421863584)

generates an array of exponential z faces 

"""
function exponential_z_faces(; Nz = 34, depth = 3000, h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
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
                        z_faces = exponential_z_faces()) 

    Nx = ceil(Int, (longitude[2] - longitude[1]) / resolution)
    Ny = ceil(Int, ( latitude[2] -  latitude[1]) / resolution)
    Nz = length(z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny, Nz),
                                            latitude,
                                            longitude,
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = z_faces)

    bathymetry(λ, φ) = neverworld_bathymetry(λ, φ, bathymetry_params; longitudinal_extent, latitude)

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry); active_cells_map = true)
end
