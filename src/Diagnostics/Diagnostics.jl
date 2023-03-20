module Diagnostics

export all_fieldtimeseries, limit_timeseries!, propagate_on_fieldtimeseries
export VolumeField, AreaField, MetricField, KineticEnergyField, time_average

using Oceananigans
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll

function propagate_on_fieldtimeseries(args...; func, nargs = 1)

    args_op   = Tuple(args[i][1] for i in 1:nargs)
    operation = func(args_op...)

    output = FieldTimeSeries{location(operation)...}(args[1].grid, args[1].times)
    
    set!(output[1], compute!(Field(operation)))

    for i in 2:length(output.times)
        @info "time $i of $(length(output.times))"
        args_op   = Tuple(args[j][i] for j in 1:nargs)
        operation = func(args_op...)

        set!(output[i], compute!(Field(operation)))
    end

    return output
end

include("spurious_mixing.jl")
include("diagnostic_fields.jl")
include("integrated_diagnostics.jl")
include("spectra.jl")

end
