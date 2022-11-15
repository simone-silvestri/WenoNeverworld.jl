using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size
using Oceananigans.Utils: instantiate


function update_simulation_clock!(simulation, init_file)
    clock = jldopen(init_file)["clock"]
    simulation.model.clock.time = clock.time	
    simulation.model.clock.iteration = clock.iteration	

    return nothing
end

function increase_simulation_Δt!(simulation; cutoff_time = 20days, new_Δt = 2minutes)
    
    function increase_Δt!(simulation)
        if simulation.model.clock.time > cutoff_time
                simulation.Δt = new_Δt
        end
    end

    simulation.callbacks[:increase_dt] = Callback(increase_Δt!, IterationInterval(1000))

    return nothing
end

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

function interpolate_per_level(old_vector, old_grid, new_grid, loc)

    H_new = halo_size(new_grid)
    H_old = halo_size(old_grid)

    Nx_new, Ny_new, _ = size(new_grid)
    Nx_old, Ny_old, _ = size(old_grid)
    Nz = size(old_vector)[3]

    k_final = Nz

    loc[2] == Face ? j_final = Ny_new + 1 : j_final = Ny_new
    
    if loc[3] == Nothing
        k_final = 1
        loc = (loc[1], loc[2], Center)
    end
    old_grid = LatitudeLongitudeGrid(CPU(), size = (Nx_old, Ny_old, 1),
                                            latitude  = (-70, 0),
                                            longitude = (-2, 62),
                                            halo = H_old,
                                            topology = (Periodic, Bounded, Bounded),
                                            z = (0, 1))


    new_grid = LatitudeLongitudeGrid(CPU(), size = (Nx_new, Ny_new, 1),
                                            latitude  = (-70, 0),
                                            longitude = (-2, 62),
                                            halo = H_new,
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

function neverworld_grid(arch, degree; H = 5)

    Nx = Int(64 / degree)
    Ny = Int(70 / degree)
    Nz = length(new_z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(arch, size = (Nx, Ny, Nz),
                                            latitude  = (-70, 0),
                                            longitude = (-2, 62),
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
