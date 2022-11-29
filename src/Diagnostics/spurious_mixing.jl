using Oceananigans.AbstractOperations: GridMetricOperation
using Oceananigans.Grids: architecture, znode
using Oceananigans.Architectures: device, device_event, arch_array

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(ρ₀ * (1 - g * b)))

function calculate_z★_diagnostics(b::FieldTimeSeries)

    times = b.times

    vol = VolumeField(b.grid)
    z★  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    total_area = sum(AreaField(b.grid))
    
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

function calculate_Γ²_diagnostics(z★::FieldTimeSeries, b::FieldTimeSeries; ρ₀ = 1000.0, g = 9.80655)
    
    times = b.times

    Γ²  = FieldTimeSeries{Center, Center, Center}(b.grid, b.times)

    for iter in 1:length(times)
        @info "time $iter of $(length(times))"

        ρ = DensityField(b[iter]; ρ₀, g)

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
    Δz       = - (z_local - z★_local) / Nint
    zrange   = z_local:Δz:z★_local

    @unroll for z in zrange
        Γ²[i, j, k] += Δz * linear_interpolate(z★_arr, ρ_arr, z)
    end
end

using Oceananigans.Fields: condition_operand

@inline function calculate_residual_MOC(v::FieldTimeSeries, b::FieldTimeSeries; blevels = 0.0:0.001:0.06)

    grid = v.grid
    arch = architecture(grid)

    Δb = blevels[2] - blevels[1]

    Nb         = length(blevels)
    Nx, Ny, Nz = size(grid)
    Nt         = length(v.times) 
    
    ψ    = zeros(Ny, Nb)
    ψint = zeros(Ny, Nb)

    for iter in Nt-1:Nt
        @info "time $iter of $(length(v.times))"
        
        for j in 1:Ny
            for i in 1:Nx, k in 1:Nz
                blev = searchsortedfirst(blevels, b[iter][i, j, k])
                ψ[j, blev] .+= v[iter][i, j, k] * Ayᶜᶠᶜ(i, j, k, v.grid) / Nt
            end
        end
    end

    ψint[:, 1] .= ψ[:, 1]

    for blev in 2:Nb
        ψint[:, blev] = ψint[:, blev-1] + Δb * ψ[:, blev]
    end

    return ψint
end

# @kernel function _calculate_ψ!(ψ, vdzdb, b, column, Nz, blevels)
#     i, j = @index(Global, NTuple)

#     if column[i, j, 1] == Nz

#         Δb = blevels[2] - blevels[1]

#         integ = Array(interior(vdzdb, i, j, :))
#         bcol  = Array(interior(b, i, j, :))

#         perm  = sortperm(bcol)
#         b_arr = bcol[perm]
#         i_arr = integ[perm]
        
#         ψ[i, j, 1] = 0.0
#         @unroll for k in 2:length(blevels)
#             ψ[i, j, k] =  ψ[i, j, k-1] - Δb * linear_interpolate(b_arr, i_arr, blevels[k])
#         end 
#     else
#         @unroll for (k, bi) in enumerate(blevels)
#             ψ[i, j, k] = 1e10
#         end 
#     end
# end

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
