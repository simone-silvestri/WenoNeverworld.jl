function integral_kinetic_energy(u::FieldTimeSeries, v::FieldTimeSeries; stride = 1, start_time = 1, end_time = length(u.times))

    energy = Float64[]
    vol = VolumeField(u.grid)

    for i in start_time:stride:end_time
        @info "integrating index $i of $end_time"
        ke = Field(u[i]^2 + v[i]^2)
        push!(energy, sum(compute!(Field(ke * vol))))
    end

    return energy
end

function ACC_transport(u::FieldTimeSeries; stride = 1, start_time = 1, end_time = length(u.times))

    transport = Float64[]
    vol = VolumeField(u.grid)

    for i in start_time:stride:end_time
        @info "integrating index $i of $end_time"
        push!(transport, sum(compute!(Field(u[i] * vol)), dims = (2, 3))[1, 1, 1])
    end

    return transport
end

function heat_content(b::FieldTimeSeries; stride = 1, start_time = 1, end_time = length(b.times))

    heat = Float64[]
    vol = VolumeField(b.grid)

    for i in start_time:stride:end_time
        @info "integrating index $i of $end_time"
        push!(heat, sum(compute!(Field(b[i] * vol))))
    end

    return heat
end

using Oceananigans.Operators: Δzᶜᶠᶜ

function calculate_eulerian_MOC(v::Field)
        
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

function calculate_eulerian_MOC(v::FieldTimeSeries)

    v̄ = time_average(v)
    ψ = calculate_MOC(v̄)
    
    return ψ
end

function time_average(field::FieldTimeSeries, iterations = 1:length(field.times))
    avg = similar(field[1])
    fill!(avg, 0)

    for t in iterations
        avg .+= field[t] ./ length(field.times)
    end

    return avg
end

function calculate_fluctuations!(fields::Dict, variables)

    for var in variables

        field_avg  = time_average(fields[var])
        field_fluc = FieldTimeSeries{location(fields[var])...}(fields[var].grid, fields[var].times)

        for t in 1:length(fields[var].times)
            fluc = compute!(Field(fields[var][t] - field_avg))
            set!(field_fluc[t], fluc)
        end

        fields[Symbol(var, :fluc)] = field_fluc
    end

    return nothing
end
