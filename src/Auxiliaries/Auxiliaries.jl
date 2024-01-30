module Auxiliaries

export cubic_interpolate, update_simulation_clock!, increase_simulation_Δt!
export parabolic_scaling, initial_buoyancy_parabola, exponential_profile
export regrid_field

using WenoNeverworld
using Oceananigans
using Oceananigans.Utils
using Oceananigans.Fields: interpolate
using Oceananigans.Grids: λnode, φnode, halo_size, on_architecture
using Oceananigans.Architectures: arch_array, architecture
using Oceananigans.Utils: instantiate
using Oceananigans.BoundaryConditions
using JLD2

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Fields: regrid!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology

include("auxiliary_functions.jl")
include("regrid_field.jl")

"""	
    function update_simulation_clock!(simulation, init_file)	

updates the `clock` of `simulation` with the time in `init_file`	
"""
function update_simulation_clock!(simulation, init_file)
    clock = jldopen(init_file)["clock"]
    simulation.model.clock.time = clock.time	
    simulation.model.clock.iteration = clock.iteration	

    return nothing
end

"""	
    function increase_simulation_Δt!(simulation; cutoff_time = 20days, new_Δt = 2minutes)

utility to update the `Δt` of a `simulation` after a certain `cutoff_time` with `new_Δt`.	
Note: this function adds a `callback` to simulation, so the order of `increase_simulation_Δt!` 	
matters (i.e. the `Δt` will be updated based on the order of `increase_simulation_Δt!` specified)	
"""	
function increase_simulation_Δt!(simulation; cutoff_time = 20days, new_Δt = 2minutes)
    
    counter = 0
    for (name, callback) in simulation.callbacks
        if occursin("increase_Δt!", string(name))
            counter = max(counter, parse(Int, string(name)[end]) + 1)
        end
    end

    increase_Δt! = Symbol(:increase_Δt!, counter)

    @eval begin
        $increase_Δt!(simulation) = simulation.Δt = $new_Δt
        callback = Callback($increase_Δt!, SpecifiedTimes(cutoff_time))
    end

    simulation.callbacks[increase_Δt!] = callback

    return nothing
end

end