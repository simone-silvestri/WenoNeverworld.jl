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

    transport = []

    for (i, time) in enumerate(u.times)
        @info "integrating time $(prettytime(time)) of $(prettytime(u.times[end]))"
        push!(transport, compute!(Field(Integral(u[i], dims = (2, 3))))[1, 1, 1])
    end

    return transport
end
