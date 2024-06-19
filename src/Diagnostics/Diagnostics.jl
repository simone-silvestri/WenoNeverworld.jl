module Diagnostics

export all_fieldtimeseries, limit_timeseries!, propagate
export VolumeField, AreaField, MetricField, KineticEnergyField, time_average, PotentialEnergyField

export KineticEnergy, VerticalVorticity, PotentialVorticity, DeformationRadius, Stratification

using Oceananigans
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils
using Oceananigans.Fields: mean
using Oceananigans.Grids: halo_size
using Oceananigans.OutputReaders: OnDisk
using JLD2
using Oceananigans.Fields: default_indices

function propagate(fields...; func)
			      
    fields_op = Tuple(retrieve_operand(field, 1) for field in fields)
    operation = func(fields_op...)

    if operation isa Number || operation isa Array 
       field_output = []
       for i in 2:length(fields[1].times)
          @info "propagate on index $i"
          fields_op = retrieve_operand.(fields, i)
          operation = func(fields_op...)
          push!(field_output, operation)
       end
       return field_output
    end

    field1 = Field(operation)
    field_output = FieldTimeSeries{location(field1)...}(fields[1].grid, fields[1].times, indices = field1.indices)

    set!(field_output, operation, 1)

    for i in 2:length(field_output.times)
        @info "propagating on index $i"
        fields_op = retrieve_operand.(fields, i)
        operation = func(fields_op...)
        set!(field_output, operation, i)
    end

    return field_output
end

import Oceananigans.Fields: set!

set!(fts::FieldTimeSeries, arr::Array, i::Int) = set!(fts[i], arr)

retrieve_operand(f::Number, i)          = f
retrieve_operand(f::Field, i)           = f
retrieve_operand(f::FieldTimeSeries, i) = f[i]

include("load_data.jl")
include("spurious_mixing.jl")
include("diagnostic_fields.jl")
include("integrated_diagnostics.jl")
include("spectra.jl")
include("compress_restart_files.jl")

end
