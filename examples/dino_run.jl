using WenoNeverworld
using WenoNeverworld.NeverworldBoundaries
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = CPU()

# The resolution in degrees
degree_resolution = 1/4

grid = NeverworldGrid(degree_resolution; arch)

# Simulation parameters
Δt        = 10minutes
stop_time = 200years

# Latitudinal wind stress acting on the zonal velocity
# a piecewise-cubic profile interpolated between
# x = φs (latitude) and y = τs (stress)
φs = (-70.0, -45.0, -15.0,  0.0,  15.0, 45.0, 70.0)
τs = (  0.0,   0.2,  -0.1, -0.01, -0.1,  0.1,  0.0)
wind_stress = WindStressBoundaryCondition(; φs, τs)

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
    T★   = ifelse(φ > 0, p.Tₙ, p.Tₛ)
    mask = sin(π * (φ + p.φₙ) / (p.φₙ - p.φₛ))

    return T★ + (p.Tₑ - T★) * mask # Need to add the time-dependent part!
end

# Idem for the salinity restoring profile
@inline function salinity_profile(λ, φ, z, t, p)
    S★   = ifelse(φ > 0, p.Sₙ, p.Sₛ)
    mask = (1 + cos(2π * φ / (p.φₙ - p.φₛ))) / 2

    return S★ + (p.Sₑ - S★) * mask - 1.25 * exp(- φ^2 / 7.5^2)
end

# The temperature flux imposed by solar radiation 
@inline function solar_flux(λ, φ, z, t, p)
    # t is in seconds, convention is that 0 is the 1st of January
    time_in_days = t / 86400
    day_of_the_year = mod(time_in_days, 365)
    solar_heat_flux = 230 * cos(π / 180 * (φ - 23.5 * cos(π * (day_of_the_year + 189) / 180)))
    return -  solar_heat_flux / p.ρ₀ / p.cₚ
end

parameters = (; Tₙ = 5.0, 
                Tₛ = - 0.5,
                Sₙ = 35.1,
                Sₛ = 35.0,
                Tₑ = 27.0,
                Sₑ = 37.25,
                φₛ = -70.0,
                φₙ = 70.0,
                ρ₀ = 1020.0,
                cₚ = 3991.0)

temperature_bc = HaneyBoundaryCondition(; restoring_profile = temperature_profile, 
                                          varname = Temperature(),
                                          flux = solar_flux,
                                          pumping_velocity = 5 / 10days,
                                          parameters)

salinity_bc = HaneyBoundaryCondition(; restoring_profile = salinity_profile, 
                                       varname = Salinity(),
                                       pumping_velocity = 5 / 10days,
                                       parameters)

tracer_boundary_conditions = (; T = temperature_bc,
                                S = salinity_bc)

# Equation of state: we use the TEOS10 equation of state
equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

# Define the initial conditions, we start with a constant salinity and
# a linear stratification in temperature going from 10 to the maximum
# temperature (27ᵒ C)
@inline initial_salinity(λ, φ, z)    = 35
@inline initial_temperature(λ, φ, z) = (grid.Lz + z) / grid.Lz * (27 - 10) + 10 # Remember! z is negative - Lz : 0

initial_conditions = (T = initial_temperature,
                      S = initial_salinity)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              wind_stress,
                                              buoyancy,
                                              tracers = (:T, :S),
                                              initial_conditions,
                                              tracer_boundary_conditions)
                                              
model = simulation.model

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

