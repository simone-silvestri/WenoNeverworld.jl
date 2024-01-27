using Oceananigans.Fields: interpolate
using Oceananigans.Grids: λnode, φnode, halo_size, on_architecture
using Oceananigans.Utils: instantiate
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations: DistributedGrid, reconstruct_global_grid

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Fields: regrid!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology

""" 
    function cubic_interpolate(x, x1, x2, y1, y2, d1, d2)

returns a cubic function between points `(x1, y1)` and `(x2, y2)` with derivative `d1` and `d2`
"""
@inline function cubic_interpolate(x; x₁, x₂, y₁, y₂, d₁ = 0, d₂ = 0)
    A = [ x₁^3 x₁^2 x₁ 1.0
          x₂^3 x₂^2 x₂ 1.0
          3*x₁^2 2*x₁ 1.0 0.0
          3*x₂^2 2*x₂ 1.0 0.0]
          
    b = [y₁, y₂, d₁, d₂]

    coeff = A \ b

    return coeff[1] * x^3 + coeff[2] * x^2 + coeff[3] * x + coeff[4]
end

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

# Disclaimer: the `_propagate_field!` implementation is copied from https://github.com/CliMA/ClimaOcean.jl/pull/60
@kernel function _propagate_field!(field, tmp_field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
        nb = (nw, ne, nn, ns)

        counter = 0
        cumsum  = 0.0 

        @unroll for n in nb
            counter += ifelse(isnan(n), 0, 1)
            cumsum  += ifelse(isnan(n), 0, n)
        end

        tmp_field[i, j, k] = ifelse(cumsum == 0, NaN, cumsum / counter)
    end
end

@kernel function _substitute_values!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    @inbounds substitute = isnan(field[i, j, k])
    @inbounds field[i, j, k] = ifelse(substitute, tmp_field[i, j, k], field[i, j, k])
end

@kernel function _nans_at_zero!(field)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(field[i, j, k] == 0, NaN, field[i, j, k])
end

propagate_horizontally!(field, ::Nothing; kw...) = nothing

""" 
    propagate_horizontally!(field; max_iter = Inf)

propagate horizontally a field with missing values at `field[i, j, k] == 0`

disclaimer:
the `propagate_horizontally!` implementation is inspired by https://github.com/CliMA/ClimaOcean.jl/pull/60
"""
function propagate_horizontally!(field; max_iter = Inf) 
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)
    
    launch!(arch, grid, :xyz, _nans_at_zero!, field)
    fill_halo_regions!(field)

    tmp_field = deepcopy(field)

    while isnan(sum(interior(field))) && iter < max_iter
        launch!(arch, grid, :xyz, _propagate_field!,   field, tmp_field)
        launch!(arch, grid, :xyz, _substitute_values!, field, tmp_field)
        iter += 1
        @debug "propagate pass $iter with sum $(sum(parent(field)))"
    end

    GC.gc()

    return nothing
end

""" 
    continue_downwards!(field)

continue downwards a field with missing values at `field[i, j, k] == 0`

the `continue_downwards!` implementation is inspired by https://github.com/CliMA/ClimaOcean.jl/pull/60
"""
function continue_downwards!(field)
    arch = architecture(field)
    grid = field.grid
    launch!(arch, grid, :xy, _continue_downwards!, field, grid)
    return nothing
end

@kernel function _continue_downwards!(field, grid)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    @unroll for k = Nz-1 : -1 : 1
        @inbounds fill_from_above = field[i, j, k] == 0
        @inbounds field[i, j, k] = ifelse(fill_from_above, field[i, j, k+1], field[i, j, k])
    end
end

function fill_missing_values!(tracer; max_iter = Inf)

    continue_downwards!(tracer)
    propagate_horizontally!(tracer; max_iter)
    
    return tracer
end

# Regrid a field in three dimensions
function three_dimensional_regrid!(a, b)
    target_grid = a.grid isa ImmersedBoundaryGrid ? a.grid.underlying_grid : a.grid
    source_grid = b.grid isa ImmersedBoundaryGrid ? b.grid.underlying_grid : b.grid 

    topo = topology(target_grid)
    arch = architecture(target_grid)
    
    yt = cpu_face_constructor_y(target_grid)
    zt = cpu_face_constructor_z(target_grid)
    Nt = size(target_grid)

    xs = cpu_face_constructor_x(source_grid)
    ys = cpu_face_constructor_y(source_grid)
    Ns = size(source_grid)

    zsize = (Ns[1], Ns[2], Nt[3])
    ysize = (Ns[1], Nt[2], Nt[3])

    # Start by regridding in z
    @debug "Regridding in z"
    zgrid   = LatitudeLongitudeGrid(arch, size = zsize, longitude = xs, latitude = ys, z = zt, topology = topo)
    field_z = Field(location(b), zgrid)
    regrid!(field_z, zgrid, source_grid, b)

    # regrid in y 
    @debug "Regridding in y"
    ygrid   = LatitudeLongitudeGrid(arch, size = ysize, longitude = xs, latitude = yt, z = zt, topology = topo)
    field_y = Field(location(b), ygrid)
    regrid!(field_y, ygrid, zgrid, field_z)

    # Finally regrid in x
    @debug "Regridding in x"
    regrid!(a, target_grid, ygrid, field_y)

    return a
end

"""	
    function regridded_field(old_vector, old_grid, new_grid, loc)	

interpolate `old_vector` (living on `loc`) from `old_grid` to `new_grid` 	
"""
function regridded_field(old_vector, old_grid, new_grid, loc)

    # Old data
    old_field = Field(loc, old_grid)
    set!(old_field, old_vector)
    
    fill_halo_regions!(old_field)
    fill_missing_values!(old_field)

    new_field = Field(loc, new_grid)
    
    return three_dimensional_regrid!(new_field, old_field)
end