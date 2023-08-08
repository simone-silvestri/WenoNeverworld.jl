using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using CairoMakie # You have to add this to your global enviroment: `] add CairoMakie`

output_dir    = joinpath(@__DIR__, "./")
@show output_prefix = output_dir * "/neverworld_quarter_resolution"

arch = GPU()

# The resolution in degrees
degree_resolution = 1/4

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

# Wanna use a different profile? Try this:
# @inline seasonal_cosine_scaling(y, t) = cos(π * y / 70) * sin(2π * t / 1year)
# buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(ΔB = 0.06, λ = 7days, func = seasonal_cosine_scaling)    

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              wind_stress,
                                              buoyancy_relaxation,
                                              interp_init,
                                              init_file)
                                              
model = simulation.model

# Let's visualize our boundary conditions!
cpu_grid = on_architecture(CPU(), grid)
φ = φnodes(cpu_grid.underlying_grid, Center(), Center(), Center())
τ_bcs = Array(model.velocities.u.boundary_conditions.top.condition.func.stress)
b_bcs = zeros(length(τ_bcs))
b = Array(interior(model.tracers.b))
for j in 1:grid.Ny
    b_bcs[j] = buoyancy_relaxation(grid.Nx÷2, j, grid, model.clock, (; b))
end

fig = Figure()
ax  = Axis(fig[1, 1], title = L"\text{Wind stress profile}")
lines!(ax, φ, - τ_bcs .* 1000, linewidth = 5) # (we re-convert the wind stress form kg/m² to Nm)

ax  = Axis(fig[1, 2], title = L"\text{Restoring buoyancy flux}")
lines!(ax, φ, - b_bcs, linewidth = 5)

CairoMakie.save("boundary_conditions.png", fig)

# Let's plot the initial conditions to make sure they are reasonable
λ = λnodes(cpu_grid.underlying_grid, Center(), Center(), Center())
z = znodes(cpu_grid.underlying_grid, Center(), Center(), Center())

fig = Figure(resolution = (800, 2000))
ax  = Axis(fig[1:4, 1], title = "Surface initial buoyancy")
hm  = heatmap!(ax, λ, φ, b[:, :, grid.Nz], colormap = :thermometer)
cb  = Colorbar(fig[1:4, 2], hm) 

ax  = Axis(fig[5, 1], title = "Mid longitude buoyancy")
hm  = heatmap!(ax, φ, z, b[grid.Nx ÷ 2, :, :], colormap = :thermometer)
ct  = contour!(ax, φ, z, b[grid.Nx ÷ 2, :, :], levels = range(0, 0.06, length = 10), color = :black)
cb  = Colorbar(fig[5, 2], hm) 

CairoMakie.save("initial_conditions.png", fig)

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
checkpoint_outputs!(simulation, output_prefix)

# initializing the time for wall_time calculation
@info "Running with Δt = $(prettytime(simulation.Δt))"
run_simulation!(simulation; interp_init, init_file)

