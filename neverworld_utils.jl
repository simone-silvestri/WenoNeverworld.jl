
using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode
using Oceananigans.Utils: instantiate

function check_zeros(grid, old_array, loc; max_passes = 10)
    Nx, Ny, Nz = size(grid)

    zc = grid.zᵃᵃᶜ[1:Nz]
    field = Field(loc, grid)
    set!(field, old_array)
    fill_halo_regions!(field)
    
    condition = false
    pass = 1
    while condition == false && pass < max_passes + 1
        @info "Pass number $pass"
        condition = true    
        for k in 1:Nz
            @info "we are at k = $k"
            for i in 1:Nx, j in 1:Ny
                if zc[k] > bathymetry[i, j]
                    if old_array[i, j, k] == 0
                        condition = false
                        @info "There is a zero! at $i, $j with z = $(zc[k]) and bat = $(bathymetry[i, j])"
                        neigbours = [field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k], field[i, j, k + 1]]
                        non_null  = Int.(neigbours .!= 0)
                        @info "Old array: $(old_array[i, j, k]), sum neigbours $(sum(neigbours))"
                        if sum(non_null) != 0 
                            old_array[i, j, k] = dot(non_null, neigbours) / sum(non_null)
                        end
                    end
                end
            end
        end
        pass += 1
    end

    return old_array
end

function interpolate_per_level(old_vector, old_degree, new_degree, loc, H, H_orig)
    
    Bλ_old = (H_orig + 2) * old_degree
    Bλ_new = (H      + 2) * new_degree

    Nx_old = Int((60 + 2Bλ_old)  / old_degree)
    Ny_old = Int(70 / old_degree)
    Nx_new = Int((60 + 2Bλ_new)  / new_degree)
    Ny_new = Int(70 / new_degree)
    Nz = 48

    loc[3] == Face ? k_final = Nz + 1 : k_final = Nz
    loc[2] == Face ? j_final = Ny_new + 1 : j_final = Ny_new
    
    if loc[3] == Nothing
        k_final = 1
        loc = (loc[1], loc[2], Center)
    end
    old_grid = LatitudeLongitudeGrid(CPU(), size = (Nx_old, Ny_old, 1),
                                            latitude  = (-70, 0),
                                            longitude = (-Bλ_old, 60+Bλ_old),
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = (0, 1))


    new_grid = LatitudeLongitudeGrid(CPU(), size = (Nx_new, Ny_new, 1),
                                            latitude  = (-70, 0),
                                            longitude = (-Bλ_new, 60+Bλ_new),
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = (0, 1))

    old_field  = Field(loc, old_grid)
    new_vector = zeros(Nx_new, j_final, k_final)

    for k in 1:k_final
        set!(old_field, old_vector[:, :, k])
        for i in 1:Nx_new, j in 1:j_final
            new_vector[i, j, k] = interpolate(old_field, xnode(loc[1](), i, new_grid), ynode(loc[2](), j, new_grid), new_grid.zᵃᵃᶜ[1])
        end
    end

    return new_vector
end

function neverworld_grid(arch, degree; H = 4)

    Nx = Int(72  / degree)
    Ny = Int(70 / degree)
    Nz = 48

    σ = 1.04 # linear stretching factor
    Δz_center_linear(k) = Lz * (σ - 1) * σ^(Nz - k) / (σ^Nz - 1) # k=1 is the bottom-most cell, k=Nz is the top cell
    linearly_spaced_faces(k) = k==1 ? -Lz : - Lz + sum(Δz_center_linear.(1:k-1))

    @show underlying_grid = LatitudeLongitudeGrid(arch, size = (Nx, Ny, Nz),
                                            latitude  = (-70, 0),
                                            longitude = (-6, 66),
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = linearly_spaced_faces)

    λ_grid = underlying_grid.λᶜᵃᵃ[1:Nx]
    φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

    bathy = zeros(Nx, Ny)
    for (i, λ) in enumerate(λ_grid), (j, φ) in enumerate(φ_grid)
        bathy[i, j] = bathymetry(λ, φ)
    end

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathy))
end

function refill_initial_condition(array, grid)
    Nx, Ny, Nz = size(grid)
    nx, ny, nz = size(array)
    new_array  = zeros(Nx, ny, nz)

    diff = Int((Nx - nx)/2)

    new_array[diff+1:diff+nx, :, :] .= array
    new_array[1:diff, :, :]         .= array[nx-diff+1:nx, :, :]
    new_array[diff+nx+1:Nx, :, :]   .= array[diff+1:2diff, :, :]

    return new_array
end