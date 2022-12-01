using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using WenoNeverworld: bathymetry_with_ridge
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, HorizontalDivergenceScalarDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalDivergenceScalarBiharmonicDiffusivity
using WenoNeverworld: geometric_νhb

output_dir    = joinpath(@__DIR__, "../files_four_isopycnal_new_bathy")
@show output_prefix = output_dir * "/neverworld_quarter_isopycnal"

arch   = GPU()
old_degree = 1
new_degree = 1/4

orig_grid = NeverworldGrid(arch, old_degree; longitude = (-5, 65)) 
grid      = NeverworldGrid(arch, new_degree)

# Extend the vertical advection scheme
interp_init = false
init_file   = "files_four_isopycnal_new_bathy/neverworld_quarter_isopycnal_checkpoint_iteration3742601.jld2" 

# Simulation parameters
Δt        = 10minutes
stop_time = 100years

using Oceananigans.Advection: compute_reconstruction_coefficients
import Oceananigans.Advection: Centered

function Centered(FT::DataType = Float64; grid = nothing, order = 2) 

    if !(grid isa Nothing) 
        FT = eltype(grid)
    end

    mod(order, 2) != 0 && throw(ArgumentError("Centered reconstruction scheme is defined only for even orders"))

    N  = Int(order ÷ 2)
    if N > 1
        coefficients = compute_reconstruction_coefficients(grid, FT, :Centered; order)
        buffer_scheme = Centered(FT; grid, order = order - 2)
    else
        coefficients    = Tuple(nothing for i in 1:6)
        buffer_scheme = nothing
    end
    return Centered{N, FT}(coefficients..., buffer_scheme)
end

upwind_advection = WENO(grid.underlying_grid)
center_advection = Centered(grid.underlying_grid, order = 6)

include("../src/isopycnally_rotated_upwinding.jl")

tracer_advection = IsopycnallyRotatedUpwindScheme(upwind_advection, center_advection)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, orig_grid, Δt, stop_time, interp_init, init_file, tracer_advection)

increase_simulation_Δt!(simulation, cutoff_time = 50days,  new_Δt = 5.0minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days, new_Δt = 7.5minutes)
increase_simulation_Δt!(simulation, cutoff_time = 300days, new_Δt = 10minutes)

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = false
standard_outputs!(simulation, output_prefix; overwrite_existing)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)
