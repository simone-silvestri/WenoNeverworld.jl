using Oceananigans.Fields: interpolate!

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
        @info "propagate pass $iter with sum $(sum(parent(field)))"
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

    topo = topology(a.grid)
    arch = architecture(a.grid)
    
    yt = cpu_face_constructor_y(a.grid)
    zt = cpu_face_constructor_z(a.grid)
    Nt = size(a.grid)

    xs = cpu_face_constructor_x(b.grid)
    ys = cpu_face_constructor_y(b.grid)
    Ns = size(b.grid)

    zsize = (Ns[1], Ns[2], Nt[3])
    ysize = (Ns[1], Nt[2], Nt[3])

    # Start by regridding in z
    @debug "Regridding in z"
    zgrid   = LatitudeLongitudeGrid(arch, size = zsize, longitude = xs, latitude = ys, z = zt, topology = topo)
    field_z = Field(location(b), zgrid)
    interpolate!(field_z, b)

    # regrid in y 
    @debug "Regridding in y"
    ygrid   = LatitudeLongitudeGrid(arch, size = ysize, longitude = xs, latitude = yt, z = zt, topology = topo)
    field_y = Field(location(b), ygrid)
    interpolate!(field_y, field_z)

    # Finally regrid in x
    @debug "Regridding in x"
    interpolate!(a, field_y)

    return a
end

"""	
    function regridded_field(old_vector, old_grid, new_grid, loc)	

interpolate `old_vector` (living on `loc`) from `old_grid` to `new_grid` 	
"""
function regrid_field(old_vector, old_grid, new_grid, loc)

    source_grid = old_grid isa ImmersedBoundaryGrid ? old_grid.underlying_grid : old_grid
    target_grid = new_grid isa ImmersedBoundaryGrid ? new_grid.underlying_grid : new_grid 

    # Old data
    old_field = Field(loc, source_grid)
    set!(old_field, old_vector)
    
    fill_halo_regions!(old_field)
    fill_missing_values!(old_field)

    new_field = Field(loc, target_grid)
    
    return three_dimensional_regrid!(new_field, old_field)
end

