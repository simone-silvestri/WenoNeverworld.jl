using WenoNeverworld
using WenoNeverworld.NeverworldGrids
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using CUDA
using LaTeXStrings

CUDA.device!(3)

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/storage4/WenoNeverworldData/quarter_degree_interp/"
@show output_prefix = output_dir * "weno_quarter_interp" 

arch = GPU()

# The resolution in degrees
degree_resolution = 1/4

z_faces = exponential_z_faces(; Nz = 35, depth = 3000)
grid = NeverworldGrid(degree_resolution; arch, z_faces)

# Do we need to interpolate? (interp_init) If `true` from which file?
interp_init = false # If interpolating from a different grid: `interp_init = true`
init_file   = "/storage4/WenoNeverworldData/half_degree_new/" * "weno_half_checkpoint_iteration42167520.jld2" # To restart from a file: `init_file = /path/to/restart`

# Simulation parameters
Δt        = 15minutes
stop_time = 2000years

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
# buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(seasonal_cosine_scaling; ΔB = 0.06, λ = 7days)    

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(grid; Δt, stop_time,
                                              wind_stress,
                                              buoyancy_relaxation,
                                              interp_init,
                                              init_file, vertical_diffusivity = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=3e-5))
                                              
model = simulation.model

# Add outputs (check other outputs to attach in `src/neverworld_outputs.jl`)
#checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = 5years)

# initializing the time for wall_time calculation
#@info "Running with Δt = $(prettytime(simulation.Δt))"
#run_simulation!(simulation; interp_init, init_file)

cpu_grid = on_architecture(CPU(), grid)
φ = φnodes(cpu_grid.underlying_grid, Center(), Center(), Center())
τ_bcs = Array(model.velocities.u.boundary_conditions.top.condition.func.stress)
#b_bcs = zeros(length(τ_bcs))
b = Array(interior(model.tracers.b))

# Let's plot the initial conditions to make sure they are reasonable
λ = λnodes(cpu_grid.underlying_grid, Center(), Center(), Center())
z = znodes(cpu_grid.underlying_grid, Center(), Center(), Center())


using CairoMakie
#using CairoMakie.AbstractPlotting.MakieLayout
#using Makie
f = Figure(resolution = (600, 800))

ax2 = Axis(f[1, 1], xticklabelcolor = :navyblue, xtickcolor = :navyblue, xlabel="Longitude [∘]", ylabel="Latitude [∘]", xticklabelsize = 20, xlabelsize=20, titlesize=20, aspect=0.5, xticks=0:20:60, xaxisposition = :bottom, yticks=-70:20:70, yticklabelsize = 20, ylabelsize=20, xlabelcolor = :navyblue)
ax1 = Axis(f[1, 1], xticklabelcolor = :red3, xlabelcolor = :red3, xtickcolor = :red3, xlabel="τ [Pa]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=-0.1:0.1:0.2, yticks=-70:20:70, xaxisposition = :top)

linkyaxes!(ax1, ax2)
hideydecorations!(ax1)
hidespines!(ax1)
ylims!(-70, 70)
#hideydecorations!(ax2)
#xaxis_top!(ax1)

hm = heatmap!(ax2, λ, φ, b[:, :, grid.Nz], colormap = :deep)
lines!(ax1, - τ_bcs .* 1000, φ, linewidth = 4, color = :red3) 
cb = CairoMakie.Colorbar(f[1, 2], hm, ticklabelsize = 20, width = 30, label = L"[m/s^2]", labelsize = 20)
# (we re-convert the wind stress form kg/m² to Nm)
f
#save("example/buoyancy_wind_plot_final.png", f)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("buoyancy_wind_plot1.png", f, px_per_unit = 3)

