using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Utils

function all_fieldtimeseries(file; variables = ("u", "v", "w", "b"))

    fields = Dict()

    for var in variables
        fields[Symbol(var)] = FieldTimeSeries(file, var; backend=OnDisk(), architecture=CPU())
    end

    return fields
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

function calculate_MOC(v::Field)
    
    dz  = MetricField((Center, Face, Center), v.grid, Oceananigans.AbstractOperations.Δz)
    moc = compute!(Field(dz * v))


    return moc
end


