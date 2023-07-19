using Oceananigans.Fields: interpolate
using Oceananigans.Grids: λnode, φnode, halo_size
using Oceananigans.Utils: instantiate
using Oceananigans.BoundaryConditions

"""	
    function update_simulation_clock!(simulation, init_file)	

updates the `clock` of `simulation` with the time in `init_file`	
"""
function update_simulation_clock!(simulation, init_file)
    clock = jldopen(init_file)["clock"]
    simulation.model.clock.time = clock.time	
    simulation.model.clock.iteration = clock.iteration	

    return nothing
end

"""	
    function increase_simulation_Δt!(simulation; cutoff_time = 20days, new_Δt = 2minutes)

utility to update the `Δt` of a `simulation` after a certain `cutoff_time` with `new_Δt`.	
Note: this function adds a `callback` to simulation, so the order of `increase_simulation_Δt!` 	
matters (i.e. the `Δt` will be updated based on the order of `increase_simulation_Δt!` specified)	
"""	
function increase_simulation_Δt!(simulation; cutoff_time = 20days, new_Δt = 2minutes)
    
    counter = 0
    for (name, callback) in simulation.callbacks
        if occursin("increase_Δt!", string(name))
            counter = max(counter, parse(Int, string(name)[end]) + 1)
        end
    end

    increase_Δt! = Symbol(:increase_Δt!, counter)

    @eval begin
        $increase_Δt!(simulation) = simulation.Δt = $new_Δt
        callback = Callback($increase_Δt!, SpecifiedTimes(cutoff_time))
    end

    simulation.callbacks[increase_Δt!] = callback

    return nothing
end

"""	
    function interpolate_per_level(old_vector, old_grid, new_grid, loc)	

interpolate `old_vector` (living on `loc`) from `old_grid` to `new_grid` 	
Note: The z-levels of `old_grid` and `new_grid` should be the same!!	
"""
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
                                            longitude = (-5, 65),
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
        fill_halo_regions!(old_field)
        for i in 1:Nx_new, j in 1:j_final
            new_vector[i, j, k] = interpolate(old_field,  λnode(i, new_grid, loc[1]()), φnode(j, new_grid, loc[2]()), new_grid.zᵃᵃᶜ[1])
        end
    end

    return new_vector
end

