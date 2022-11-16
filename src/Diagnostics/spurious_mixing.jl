using Oceananigans.AbstractOperations: GridMetricOperation

VolumeField(loc, grid) = compute!(Field(GridMetricOperation(loc, Oceananigans.AbstractOperations.volume, grid)))
AreaField(loc, grid)   = Compute!(Field(GridMetricOperation(loc, Oceananigans.AbstractOperations.Az, grid)))

MetricField(loc, grid, metric)   = Compute!(Field(GridMetricOperation(loc, metric, grid)))

function calculate_z★_diagnostics(b::FieldTimeSeries)

    times = b.times

    vol = VolumeField((Center, Center, Center), b.grid)
    z★  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    total_area = sum(AreaField((Center, Center, Nothing), b.grid))
    
    for iter in 1:length(times)
       @info "time $iter of $(length(times))"
       calculate_z★!(z★[iter], b[iter], vol, total_area)
    end
        
    return z★
end

function calculate_z★!(z★::Field, b::Field, vol, total_area)
    grid = b.grid
    arch = architecture(grid)

    b_arr = Array(interior(b))[:]
    v_arr = Array(interior(vol))[:]

    perm           = sortperm(b_arr)
    sorted_b_field = b_arr[perm]
    sorted_v_field = v_arr[perm]
    integrated_v   = cumsum(sorted_v_field)    

    z★_event = launch!(arch, grid, :xyz, _calculate_z★, z★, b, sorted_b_field, integrated_v; dependencies = device_event(arch))
    wait(device(arch), z★_event)

    z★ ./= total_area

    return nothing
end

@kernel function _calculate_z★(z★, b, b_sorted, integrated_v)
    i, j, k = @index(Global, NTuple)
    bl  = b[i, j, k]
    i₁  = searchsortedfirst(b_sorted, bl)
    z★[i, j, k] = integrated_v[i₁] 
end

function calculate_Γ²_diagnostics(z★::FieldTimeSeries, b::FieldTimeSeries)
    
    times = b.times

    Γ²  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    for iter in 1:length(times)
        @info "time $iter of $(length(times))"
        calculate_Γ²!(Γ²[iter], z★[iter], b[iter])
    end
         
    return Γ²
end

function calculate_Γ²!(Γ², z★, b)
    perm   = sortperm(Array(interior(z★))[:])

    b_arr  = (Array(interior(b))[:])[perm]
    z★_arr = (Array(interior(z★))[:])[perm]

    Γ²_event = launch!(arch, grid, :xyz, _calculate_Γ², Γ², z★, z★_arr, b_arr, grid; dependencies = device_event(arch))
    wait(device(arch), Γ²_event)

    return nothing
end

@kernel function _calculate_Γ²(Γ², z★, z★_arr, b_arr, grid)
    i, j, k = @index(Global, NTuple)

    Nint = 10.0
     
    Γ²[i, j, k] = 0.0
         
    z_local  = znode(Center(), k, grid) + grid.Lz
    z★_local = z★[i, j, k] 
    Δz       = (z_local - z★_local) / Nint
    zrange   = z★_local:Δz:z_local

    @unroll for z in zrange
        Γ²[i, j, k] += Δz * linear_interpolate(z★_arr, b_arr, z)
    end
end

@inline function linear_interpolate(x, y, x₀)
    i₁ = searchsortedfirst(x, x₀)
    i₂ =  searchsortedlast(x, x₀)

    @inbounds y₂ = y[i₂]
    @inbounds y₁ = y[i₁]

    @inbounds x₂ = x[i₂]
    @inbounds x₁ = x[i₁]

    if x₁ == x₂
        return y₁
    else
        return (y₂ - y₁) / (x₂ - x₁) * (x₀ - x₁) + y₁
    end
end