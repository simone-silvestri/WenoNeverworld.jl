module Diagnostics

using Oceananigans
using KernelAbstractions: @kernel, @index 
using KernelAbstractions.Extras.LoopInfo: @unroll

include("spurious_mixing.jl")
include("integrated_diagnostics.jl")
include("spectra.jl")

end