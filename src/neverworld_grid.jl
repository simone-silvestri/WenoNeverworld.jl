using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size

z_faces = [-4000.0, -3740.0, -3422.0, -3126.0, -2854.0, -2604.0, -2378.0, -2174.0, -1993.0, -1834.0, -1695.0, -1572.0,
 -1461.0, -1356.0, -1255.0, -1155.0, -1056.0, -958.0, -861.0, -767.0, -677.0, -592.0, -513.0, -441.0, -378.0, -323.0, -276.0, -238.0,
 -207.0, -182.0, -162.0, -146.0, -133.0, -121.0, -110.0, -100.0, -90.0, -80.0, -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, -10.0, 0.0]

all_z_faces = zeros(length(z_faces)*2-1) 
for (i, z) in enumerate(all_z_faces)
    if mod(i, 2) == 0
        all_z_faces[i] = 0.5 * (z_faces[i ÷ 2] + z_faces[i ÷ 2 + 1])
    else
        all_z_faces[i] = z_faces[Int(i ÷ 2) + 1]
    end
end

new_z_faces = vcat(z_faces[1:22], all_z_faces[44:end])

using Oceananigans.ImmersedBoundaries: PartialCellBottom

function NeverworldGrid(arch, degree; H = 5, longitude = (-2, 62), latitude = (-70, 0), bathymetry = bathymetry_without_ridge)

    Nx = Int((longitude[2] - longitude[1]) / degree)
    Ny = Int((latitude[2]  - latitude[1]) / degree)
    Nz = length(new_z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz),
                                            latitude,
                                            longitude,
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = new_z_faces)

    λ_grid = underlying_grid.λᶜᵃᵃ[1:Nx]
    φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

    bathy = zeros(Nx, Ny)
    for (i, λ) in enumerate(λ_grid), (j, φ) in enumerate(φ_grid)
        bathy[i, j] = bathymetry(λ, φ)
    end

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathy))
end
