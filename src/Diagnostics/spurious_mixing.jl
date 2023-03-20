using Oceananigans.AbstractOperations: GridMetricOperation
using Oceananigans.Grids: architecture, znode
using Oceananigans.Architectures: device, arch_array

MetricField(loc, grid, metric; indices = default_indices(3)) = compute!(Field(GridMetricOperation(loc, metric, grid); indices))

VolumeField(grid, loc=(Center, Center, Center);  indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume; indices)
  AreaField(grid, loc=(Center, Center, Nothing); indices = default_indices(3)) = MetricField(loc, grid, Oceananigans.AbstractOperations.Az; indices)

DensityField(b::Field; ρ₀ = 1000.0, g = 9.80655) = compute!(Field(ρ₀ * (1 - g * b)))

function HeightField(grid, loc = (Center, Center, Center))  

    zf = Field(loc, grid)

    for k in 1:size(zf, 3)
        interior(zf, :, :, k) .= znode(loc[3](), k, grid)
    end

    return zf
end

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

    launch!(arch, grid, :xyz, _calculate_z★, z★, b, sorted_b_field, integrated_v)
    
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

    launch!(arch, grid, :xyz, _calculate_Γ², Γ², z★, z★_arr, ρ_arr, grid)

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

@inline function calculate_κeff(b::FieldTimeSeries, ψint; blevels = collect(0.0:0.001:0.06))

    grid = b.grid
    arch = architecture(grid)

    Nb         = length(blevels)
    Nx, Ny, Nz = size(grid)
    Nt         = length(b.times) 
    
    strat       = [zeros(Nx, Ny, Nb) for iter in 1:Nt]
    stratint    = [zeros(Nx, Ny, Nb) for iter in 1:Nt]
    stratavgint = zeros(Ny, Nb)

    ψint2 = zeros(Ny, Nb)

    for iter in 1:Nt
        @info "time $iter of $(length(b.times))"

        bz = compute!(Field(∂z(b[iter])))
        launch!(arch, grid, :xy, _cumulate_stratification!, stratint[iter], strat[iter], bz, b[iter], blevels, grid, Nz)
    end

    for iter in 1:Nt
        for i in 20:220
            stratavgint .+= stratint[iter][i, :, :] / Nt
        end
    end

    Δb = blevels[2] - blevels[1]

    for j in 1:Ny
        ψint2[j, 1] = Δb * ψint[j, 1]
        for blev in 2:Nb
            ψint2[j, blev] = ψint2[j, blev-1] + Δb * ψint[j, blev]
        end
    end

    return ψint2 ./ stratavgint
end

@kernel function _cumulate_stratification!(stratint, strat, bz, b, blevels, grid, Nz)
    i, j = @index(Global, NTuple)

    Nb = length(blevels)
    Δb = blevels[2] - blevels[1]

    @unroll for k in 1:Nz
        if b[i, j, k] < blevels[end]
            blev = searchsortedfirst(blevels, b[i, j, k])
            strat[i, j, blev] += bz[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) 
        end
    end

    stratint[i, j, 1] = Δb * strat[i, j, 1]
    bmax = maximum(b[i, j, :])
    @unroll for blev in 2:Nb
        if bmax > blevels[blev]
            stratint[i, j, blev] = stratint[i, j, blev-1] + Δb * strat[i, j, blev]
        end
    end
end

@inline function calculate_residual_MOC(v::FieldTimeSeries, b::FieldTimeSeries; blevels = collect(0.0:0.001:0.06))

    grid = v.grid
    arch = architecture(grid)

    Nb         = length(blevels)
    Nx, Ny, Nz = size(grid)
    Nt         = length(v.times) 
    
    ψ       = [zeros(Nx, Ny, Nb) for iter in 1:Nt]
    ψint    = [zeros(Nx, Ny, Nb) for iter in 1:Nt]
    ψavgint = zeros(Ny, Nb)

    for iter in 1:Nt
        @info "time $iter of $(length(v.times))"
        launch!(arch, grid, :xy, _cumulate_v_velocities!, ψint[iter], ψ[iter], b[iter], v[iter], blevels, grid, Nz)
    end

    for iter in 1:Nt
        for i in 20:220
            ψavgint .+= ψint[iter][i, :, :] / Nt
        end
    end

    return ψavgint
end

@kernel function _cumulate_v_velocities!(ψint, ψ, b, v, blevels, grid, Nz)
    i, j = @index(Global, NTuple)

    Nb = length(blevels)
    Δb = blevels[2] - blevels[1]

    @unroll for k in 1:Nz
        if b[i, j, k] < blevels[end]
            blev = searchsortedfirst(blevels, b[i, j, k])
            ψ[i, j, blev] += v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) 
        end
    end

    ψint[i, j, 1] = Δb * ψ[i, j, 1]
    bmax = maximum(b[i, j, :])
    @unroll for blev in 2:Nb
        if bmax > blevels[blev]
            ψint[i, j, blev] = ψint[i, j, blev-1] + Δb * ψ[i, j, blev]
        end
    end

end

@inline function linear_interpolate(x, y, x₀)
    i₁ = searchsortedfirst(x, x₀)
    i₂ =  searchsortedlast(x, x₀)

    @inbounds y₂ = y[i₂]
    @inbounds y₁ = y[i₁]

    @inbounds x₂ = x[i₂]
    @inbounds x₁ = x[i₁]

    if i₁ > length(x)
        return y₂
    elseif i₁ == i₂
        isnan(y₁) && @show i₁, i₂, x₁, x₂, y₁, y₂
        return 
    else
        if isnan(y₁) || isnan(y₂) || isnan(x₁) || isnan(x₂) 
            @show i₁, i₂, x₁, x₂, y₁, y₂
        end
        return (y₂ - y₁) / (x₂ - x₁) * (x₀ - x₁) + y₁
    end
end

