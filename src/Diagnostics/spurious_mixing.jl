using Oceananigans.AbstractOperations: GridMetricOperation
using Oceananigans.Grids: architecture, znode
using Oceananigans.Architectures: device, device_event, arch_array

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(ρ₀ * (1 - g * b)))

function calculate_z★_diagnostics(b::FieldTimeSeries; ρ₀ = 1000.0, g = 9.80655)

    times = b.times

    vol = VolumeField(b.grid)
    z★  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    total_area = sum(AreaField(b.grid))
    
    for iter in 1:length(times)
       @info "time $iter of $(length(times))"

       ρ = DensityField(b; ρ₀, g)

       calculate_z★!(z★[iter], ρ, vol, total_area)
    end
        
    return z★
end

function calculate_z★!(z★::Field, ρ::Field, vol, total_area)
    grid = b.grid
    arch = architecture(grid)

    ρ_arr = Array(interior(ρ))[:]
    v_arr = Array(interior(vol))[:]

    perm           = sortperm(ρ_arr)
    sorted_ρ_field = ρ_arr[perm]
    sorted_v_field = v_arr[perm]
    integrated_v   = cumsum(sorted_v_field)    

    z★_event = launch!(arch, grid, :xyz, _calculate_z★, z★, ρ, sorted_ρ_field, integrated_v; dependencies = device_event(arch))
    wait(device(arch), z★_event)

    z★ ./= total_area

    return nothing
end

@kernel function _calculate_z★(z★, ρ, ρ_sorted, integrated_v)
    i, j, k = @index(Global, NTuple)
    ρl  = ρ[i, j, k]
    i₁  = searchsortedfirst(ρ_sorted, ρl)
    z★[i, j, k] = integrated_v[i₁] 
end

function calculate_Γ²_diagnostics(z★::FieldTimeSeries, b::FieldTimeSeries; ρ₀ = 1000.0, g = 9.80655)
    
    times = b.times

    Γ²  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    for iter in 1:length(times)
        @info "time $iter of $(length(times))"

        ρ = DensityField(b; ρ₀, g)

        calculate_Γ²!(Γ²[iter], z★[iter], ρ)
    end
         
    return Γ²
end

function calculate_Γ²!(Γ², z★, ρ)
    grid = ρ.grid
    arch = architecture(grid)

    perm   = sortperm(Array(interior(z★))[:])

    ρ_arr  = (Array(interior(ρ))[:])[perm]
    z★_arr = (Array(interior(z★))[:])[perm]

    Γ²_event = launch!(arch, grid, :xyz, _calculate_Γ², Γ², z★, z★_arr, ρ_arr, grid; dependencies = device_event(arch))
    wait(device(arch), Γ²_event)

    return nothing
end

@kernel function _calculate_Γ²(Γ², z★, z★_arr, ρ_arr, grid)
    i, j, k = @index(Global, NTuple)

    Nint = 10.0
     
    Γ²[i, j, k] = 0.0
         
    z_local  = znode(Center(), k, grid) + grid.Lz
    z★_local = z★[i, j, k] 
    Δz       = (z_local - z★_local) / Nint
    zrange   = z★_local:Δz:z_local

    @unroll for z in zrange
        Γ²[i, j, k] += Δz * linear_interpolate(z★_arr, ρ_arr, z)
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