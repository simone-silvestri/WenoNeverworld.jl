using Oceananigans.Utils
using Oceananigans.Fields: mean
using Oceananigans.Grids: halo_size
using Oceananigans.OutputReaders: OnDisk
using JLD2
using Oceananigans.Fields: default_indices

function checkpoint_fields(file)
    file = jldopen(file)
    grid = file["grid"]

    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    b = CenterField(grid)

    Hx, Hy, Hz = halo_size(grid)
    for (var, name) in zip((u, v, w, b), ("u", "v", "w", "b"))
        set!(var, file[name * "/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz])
        fill_halo_regions!(var)
    end

    return (; u, v, w, b)
end

function all_fieldtimeseries(file; variables = ("u", "v", "w", "b"))

    fields = Dict()

    for var in variables
        fields[Symbol(var)] = FieldTimeSeries(file, var; backend=OnDisk(), architecture=CPU())
    end

    return fields
end

function limit_timeseries!(fields::Dict, times)
    new_fields = Dict()

    for (key, field) in fields
        new_fields[key] = limit_timeseries!(field, times)
    end

    return new_fields
end

function limit_timeseries!(field::FieldTimeSeries, times)

    loc = location(field)
    new_field = FieldTimeSeries{loc...}(field.grid, times)

    for (idx, time) in enumerate(field.times)
        id2 = findfirst(isequal(time), times)
        if !isnothing(id2)
            set!(new_field[id2], field[idx])
        end
    end

    return new_field
end

function add_kinetic_energy_and_vorticity_timeseries!(fields::Dict)

    ζ =     FieldTimeSeries{Face, Face, Center}(fields[:u].grid, fields[:u].times)
    E = FieldTimeSeries{Center, Center, Center}(fields[:u].grid, fields[:u].times)

    for t in 1:length(E.times)
        set!(ζ[t], VerticalVorticityField(fields, t))
        set!(E[t], KineticEnergyField(fields, t))
    end

    fields[:E] = E
    fields[:ζ] = ζ

    return nothing
end

function kinetic_energy(u::FieldTimeSeries, v::FieldTimeSeries)

    energy = Float64[]

    for (i, time) in enumerate(u.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(u.times[end]))"
        ke = Field(u[i]^2 + v[i]^2)
        push!(energy, compute!(Field(Integral(ke)))[1, 1, 1])
    end

    return energy
end

function ACC_transport(u::FieldTimeSeries)

    transport = Float64[]
    vol = VolumeField(u.grid)

    for (i, time) in enumerate(u.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(u.times[end]))"
        push!(transport, sum(compute!(Field(u[i] * vol)), dims = (2, 3))[1, 1, 1])
    end

    return transport
end

function heat_content(b::FieldTimeSeries)

    heat = Float64[]
    vol = VolumeField(b.grid)

    for (i, time) in enumerate(b.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(b.times[end]))"
        push!(heat, sum(compute!(Field(b[i] * vol))))
    end

    return heat
end

using Oceananigans.Operators: Δzᶜᶠᶜ

function calculate_MOC(v::Field)
        
    v̄ = compute!(Field(Integral(v, dims = 1)))

    ψ = Field((Nothing, Face, Face), v.grid)

    for k in 2:v.grid.Nz
        dz =  Δzᶜᶠᶜ(1, 1, k-1, v.grid)
        for j in 1:size(v.grid, 2)
            ψ[1, j, k] = ψ[1, j, k - 1] + v̄[1, j, k - 1] * dz
        end
    end

    return ψ
end

function calculate_MOC(v::FieldTimeSeries)

    v̄ = time_average(v)
    ψ = calculate_MOC(v̄)
    
    return ψ
end

function time_average(field::FieldTimeSeries)
    avg = similar(field[1])
    fill!(avg, 0)

    for t in 1:length(field.times)
        avg .+= field[t] ./ length(field.times)
    end

    return avg
end