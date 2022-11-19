using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Utils
using Oceananigans.Fields: mean

function all_fieldtimeseries(file; variables = ("u", "v", "w", "b"))

    fields = Dict()

    for var in variables
        fields[Symbol(var)] = FieldTimeSeries(file, var; backend=OnDisk(), architecture=CPU())
    end

    return fields
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

    for (i, time) in enumerate(u.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(u.times[end]))"
        push!(transport, compute!(Field(Integral(u[i], dims = (2, 3))))[1, 1, 1])
    end

    return transport
end

function heat_content(b::FieldTimeSeries)

    heat = Float64[]
    
    for (i, time) in enumerate(b.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(b.times[end]))"
        push!(heat, compute!(Field(Integral(b[i])))[1, 1, 1])
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

function calculate_MOC(v::FieldTimeSeries, indices)
        
    v̄ = similar(v[1])

    fill!(v̄, 0.0)

    for idx in indices
        v̄ .+= compute!(Field(Integral(v[idx], dims = 1))) / length(indices)
    end

    ψ = Field((Nothing, Face, Face), v.grid)

    for k in 2:v.grid.Nz
        dz =  Δzᶜᶠᶜ(1, 1, k-1, v.grid)
        for j in 1:size(v.grid, 2)
            ψ[1, j, k] = ψ[1, j, k - 1] + v̄[1, j, k - 1] * dz
        end
    end

    return ψ
end


