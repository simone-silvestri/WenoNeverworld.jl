using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes
using CairoMakie

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = CPU()

# The resolution in degrees
degree_resolution = 1

grid = NeverworldGrid(degree_resolution; arch)

# Do we need to interpolate? (interp_init) If `true` from which file?
interp_init = false # If interpolating from a different grid: `interp_init = true`
init_file   = nothing # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 10minutes
stop_time = 200years

# Latitudinal wind stress acting on the zonal velocity
# a piecewise-cubic profile interpolated between
# x = φs (latitude) and y = τs (stress)
φs = (-70.0, -45.0, -15.0,  0.0,  15.0, 45.0, 70.0)
τs = (  0.0,   0.2,  -0.1, -0.02, -0.1,  0.1,  0.0)
wind_stress = WindStressBoundaryCondition(; φs, τs)

# Buoyancy relaxation profile:
# a parabolic profile between 0, at the poles, and ΔB = 0.06 at the equator
# the restoring time is λ = 7days
buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(ΔB = 0.06, λ = 7days)

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              wind_stress,
                                              buoyancy_relaxation)
                                              
model = simulation.model

# Let's visualize our boundary conditions!
φ = φnodes(grid.underlying_grid, Center(), Center(), Center())
τ = model.velocities.u.boundary_conditions.top.condition.func.stress
b = zeros(length(τ))
for j in 1:grid.Ny
    b[j] = buoyancy_relaxation(grid.Nx÷2, j, grid, model.clock, model.tracers)
end

fig = Figure()
ax  = Axis(fig[1, 1], title = L"\text{Wind stress profile}")
lines!(ax, φ, - τ .* 1000, linewidth = 5) # (we convert the wind stress in Nm)

ax  = Axis(fig[1, 2], title = L"\text{Restoring buoyancy flux}")
lines!(ax, φ, b, linewidth = 5)

CairoMakie.save("boundary_conditions.png", fig)

# Let's plot the initial conditions to make sure they are reasonable
λ = λnodes(grid.underlying_grid, Center(), Center(), Center())
z = znodes(grid.underlying_grid, Center(), Center(), Center())

b = model.tracers.b

fig = Figure()
ax  = Axis(fig[1, 1], title = "Surface initial buoyancy")
hm  = heatmap!(ax, λ, φ, Array(interior(b, :, :, grid.Nz)), colormap = :thermometer)
cb  = Colorbar(fig[1, 2], hm) 

ax  = Axis(fig[2, 1], title = "Mid longitude buoyancy")
hm  = heatmap!(ax, φ, z, Array(interior(b, grid.Nx ÷ 2, :, :)), colormap = :thermometer)
ct  = contour!(ax, φ, z, Array(interior(b, grid.Nx ÷ 2, :, :)), levels = range(0, 0.06, length = 10), color = :black)
cb  = Colorbar(fig[2, 2], hm) 

CairoMakie.save("initial_conditions.png", fig)

# Add outputs
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

