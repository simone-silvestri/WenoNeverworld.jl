using WenoNeverworld
using WenoNeverworld.NeverworldBoundaries
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using SeawaterPolynomials
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

# Equation of state: we use the TEOS10 equation of state
equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)
ρTOES10  = equation_of_state.reference_density
cTEOS10  = SeawaterPolynomials.TEOS10.teos10_reference_heat_capacity

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
    T★   = ifelse(φ > 0, p.Tⁿ, p.Tˢ)
    mask = sin(π * (φ + p.φⁿ) / (p.φⁿ - p.φˢ))

    return T★ + (p.Tᵉ - T★) * mask # Need to add the time-dependent part!
end

# Idem for the salinity restoring profile
@inline function salinity_profile(λ, φ, z, t, p)
    S★   = ifelse(φ > 0, p.Sⁿ, p.Sˢ)
    mask = (1 + cos(2π * φ / (p.φⁿ - p.φˢ))) / 2

    return S★ + (p.Sᵉ - S★) * mask - 1.25 * exp(- φ^2 / 7.5^2)
end

# The temperature flux imposed by solar radiation 
@inline function solar_flux(λ, φ, z, t, p)
    # t is in seconds, convention is that 0 is the 1st of January
    time_in_days = t / 86400
    day_of_the_year = mod(time_in_days, 365)
    solar_heat_flux = p.Q⁰ * cos(π / 180 * (φ - p.δ * cos(π * (day_of_the_year + 189) / 180)))
    return - solar_heat_flux / p.ρ⁰ / p.cᵖ
end

# Parameters to use in the boundary conditions
parameters = (; Tⁿ = 5.0,     # temperature restoring at northern boundary
                Tˢ = - 0.5,   # temperature restoring at southern boundary
                Tᵉ = 27.0,    # temperature restoring at equator
                Sⁿ = 35.1,    # salinity restoring at northern boundary 
                Sˢ = 35.0,    # salinity restoring at southern boundary
                Sᵉ = 37.25,   # salinity restoring at equator
                φˢ = -70.0,   # southnmost edge
                φⁿ = 70.0,    # northernmost edge
                ρ⁰ = ρTEOS10, # reference density
                cᵖ = cTEOS10, # reference heat capacit
                Q⁰ = 230.0,   # reference solar flux
                δ  = 23.5)    # declination angle

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

# Define the initial conditions, we start with a constant salinity and
# a linear stratification in temperature going from 10 to the maximum
# temperature (27ᵒ C). 
# TODO: David, if you have already evolved initial conditions, please use them
# so we don't have to run for 2000 years to equilibrate! You can provide them as
# arrays of size (Nx, Ny, Nz) (the size of the grid)
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
                                              
# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

