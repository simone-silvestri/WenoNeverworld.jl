using WenoNeverworld
using WenoNeverworld.Auxiliaries
using WenoNeverworld.NeverworldBoundaries
using WenoNeverworld.NeverworldGrids: dino_parameters
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnode, λnode, znode, node
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: 
                                    CATKEMixingLength, 
                                    CATKEVerticalDiffusivity
using Oceananigans.Operators
using SeawaterPolynomials
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using NCDatasets

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = GPU()

using CUDA
CUDA.device!(1)

# The resolution in degrees
resolution = 1/4
H = 10 # this should be `max(padding, 7)` (because of advection stencil) 
fill_land_in_halos = true
grid = NeverworldGrid(resolution; arch, H, fill_land_in_halos, dino_parameters(resolution)...)

# Simulation parameters (we start with 1 minute timestep and
# increase it as the simulation equilibrates)
starting_Δt = 2minutes
stop_time   = 200years

# Equation of state: we use the TEOS10 equation of state
equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)
ρTEOS10  = equation_of_state.reference_density
cTEOS10  = SeawaterPolynomials.TEOS10.teos10_reference_heat_capacity

# Boundary conditions for Salinity and Temperature, a mix of a flux and a restoring,
# also called Haney boundary conditions...

# Restoring profile for temperature as a function of λ, φ, z, t, p
# where 
# λ - longitude
# φ - latitude
# z - depth
# t - time (in seconds)
# p - additional parameters we pass to the BC constructor (a named tuple)
@inline function temperature_profile(λ, φ, z, t, p)

    time_in_days = t / 86400
    day_of_the_year = mod(time_in_days, 360)

    T★ⁿ = p.Tⁿ + p.Tₛⁿ * seasonal_cycle(day_of_the_year, 30)
    T★ˢ = p.Tˢ + p.Tₛˢ * seasonal_cycle(day_of_the_year, 30)

    T★   = ifelse(φ > 0, T★ⁿ, T★ˢ)
    mask = sin(π * (φ + p.φⁿ) / (p.φⁿ - p.φˢ))

    return T★ + (p.Tᵉ - T★) * mask 
end

@inline function seasonal_cycle(day_of_the_year, lag)
    time_max =  5 * 30 + 21 + lag  # 21th june     at 24h in hours
    time_min = 11 * 30 + 21 + lag  # 21th december        in hours

    seasonal_cycle = cos(π * (day_of_the_year - time_max) / (time_min - time_max) )
    
    return seasonal_cycle
end

# Idem for the salinity restoring profile
@inline function salinity_profile(λ, φ, z, t, p)
    S★   = ifelse(φ > 0, p.Sⁿ, p.Sˢ)
    mask = (1 + cos(2π * φ / (p.φⁿ - p.φˢ))) / 2

    return S★ + (p.Sᵉ - S★) * mask - 1.25 * exp(- φ^2 / 7.5^2)
end

@inline function solar_flux(λ, φ, z, t, p)
    # t is in seconds, convention is that 0 is the 1st of January
    time_in_days = t / 86400
    day_of_the_year = mod(time_in_days, 360)

    solar_heat_flux = p.Q⁰ * cos(π / 180 * (φ - p.δ * cos(π * (day_of_the_year + 189) / 180)))

    return solar_heat_flux * (p.ρ⁰⁻¹ * p.cᵖ⁻¹)
end

# Parameters to use in the boundary conditions
parameters = (; Tⁿ    = 5.0,         # temperature restoring at northern boundary
                Tˢ    = - 0.5,       # temperature restoring at southern boundary
                Tᵉ    = 27.0,        # temperature restoring at equator
                Tₛⁿ   = 3.0,         # seasonal temperature restoring correction north
                Tₛˢ   = 0.5,         # seasonal temperature restoring correction north
                Sⁿ    = 35.1,        # salinity restoring at northern boundary 
                Sˢ    = 35.0,        # salinity restoring at southern boundary
                Sᵉ    = 37.25,       # salinity restoring at equator
                φˢ    = -70.0,       # southnmost edge
                φⁿ    = 70.0,        # northernmost edge
                ρ⁰⁻¹  = 1 / ρTEOS10, # reciprocal reference density 
                cᵖ⁻¹  = 1 / cTEOS10, # reciprocal reference heat capacity
                Q⁰    = 230.0,       # reference solar flux
                δ     = 23.5,        # declination angle
                Rᴿ    = 0.58,        # fraction of red light
                ξᴿ    = 1 / 0.2,     # extintion length of red light m⁻¹
                ξᴮ    = 1 / 25)      # extintion length of blue light m⁻¹

temperature_bc = HaneyBoundaryCondition(; restoring_profile = temperature_profile,
                                          flux = solar_flux,
                                          varname = Temperature(),
                                          pumping_velocity = 5 / 10days,
                                          parameters)

salinity_bc = HaneyBoundaryCondition(; restoring_profile = salinity_profile, 
                                       varname = Salinity(),
                                       pumping_velocity = 5 / 10days,
                                       parameters)

tracer_boundary_conditions = (; T = temperature_bc,
                                S = salinity_bc)

# The heating imposed by penetrative solar radiation 
@inline function solar_heating(i, j, k, grid, clock, fields, p)
    λ, φ, z = node(i, j, k, grid.underlying_grid, Center(), Center(), Center())
    
    S  = solar_flux(λ, φ, z, clock.time, p)

    Sᴿ = S * p.Rᴿ
    Sᴮ = S * (1 - p.Rᴿ)

    z⁺ = znode(k+1, grid.underlying_grid, Face())
    z⁻ = znode(k,   grid.underlying_grid, Face())

    # z⁺ and z⁻ are negative nodes!
    S⁺ = Sᴿ * exp(z⁺ * p.ξᴿ) + Sᴮ * exp(z⁺ * p.ξᴮ)
    S⁻ = Sᴿ * exp(z⁻ * p.ξᴿ) + Sᴮ * exp(z⁻ * p.ξᴮ)

    return  (S⁺ - S⁻) / Δzᶜᶜᶜ(i, j, k, grid) 
end

solar_forcing = Forcing(solar_heating; discrete_form = true, parameters)

# Initial conditions
initial_conditions_data = Dataset(joinpath(@__DIR__, "TS_init_1_4degree.nc"))

initial_salinity    = reverse(PermutedDimsArray(initial_conditions_data["soce"][:, :, :], (3, 2, 1)), dims = 3) |> Array{Float32}
initial_temperature = reverse(PermutedDimsArray(initial_conditions_data["toce"][:, :, :], (3, 2, 1)), dims = 3) |> Array{Float32}

initial_conditions = (T = initial_temperature,
                      S = initial_salinity)
                      
# Add parameterizations
# horizontal_closure = NNbackscatteringClosure(; architecture = arch, weight_path = "....")
# mixing_length = CATKEMixingLength(; Cᵇ = 0.01)
# vertical_diffusivity = CATKEVerticalDiffusivity(; mixing_length)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt = starting_Δt, stop_time,
                                              buoyancy,
                                              tracers = (:T, :S),
                                            #   forcing = (; T = solar_forcing),
                                              initial_conditions,
                                              tracer_boundary_conditions)
                                 

include("propagate_initial_conditions.jl")

propagate_horizontally!(simulation.model.tracers.S)
propagate_horizontally!(simulation.model.tracers.T)

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
reduced_outputs!(simulation, output_prefix;
                 checkpoint_time = 100days,
                   snapshot_time = 1000days,
                    surface_time = 10days)

# Initialize with a small time step and increase it after the 
# inital spin up has completed
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 15minutes, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run!(simulation)

