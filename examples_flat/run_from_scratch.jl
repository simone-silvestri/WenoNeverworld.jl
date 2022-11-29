using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: VelocityStencil
using WenoNeverworld: reduced_outputs!

output_dir    = joinpath(@__DIR__, "../files_four_centered_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter_centered"

@inline function zonal_reduced_wind_stress(y, mid_wind)
    if y < -45
        return cubic_profile(y, -70.0, -45.0, 0.0, 0.2, 0.0, 0.0)
    elseif y < -30
        return cubic_profile(y, -45.0, -30, 0.2, 0.0, 0.0, 0.0)
    else
        return cubic_profile(y, 45.0, 70.0, 0.1, 0.0, 0.0, 0.0)
    end
end

interp_init = false
init_file   = nothing #"files_four_centered_new_bathy/neverworld_quarter_centered_checkpoint_iteration3742601.jld2" 

ibg = NeverWorldGrid(GPU(), 1/4, Float32; latitude = (-70, -30), longitude = (-1, 31), longitudinal_extent = 30)

Δt        = 5minutes
stop_time = 1years

tracer_advection   = WENO(ibg.underlying_grid)
momentum_advection = WENO(Float32; vector_invariant = VelocityStencil())

simulation = weno_neverworld_simulation(; grid = ibg, Δt, interp_init, init_file, stop_time, wind_stress = zonal_reduced_wind_stress, tracer_advection, momentum_advection)

increase_simulation_Δt!(simulation, cutoff_time = 30days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 60days, new_Δt = 10minutes)

reduced_outputs!(simulation, output_prefix)

run_simulation!(simulation; interp_init, init_file)