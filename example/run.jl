using WenoNeverworld
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnodes, λnodes, znodes, on_architecture
using GLMakie

using LaTeXStrings
#using CairoMakie # You have to add this to your global enviroment: `] add CairoMakie`

output_dir    = joinpath(@__DIR__, "./")
output_dir = ""
@show output_prefix = output_dir * "" #"WenoNeverworldData/weno_fourth"

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
# buoyancy_relaxation = BuoyancyRelaxationBoundaryCondition(seasonal_cosine_scaling; ΔB = 0.06, λ = 7days)    

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
#b_bcs = zeros(length(τ_bcs))
b = Array(interior(model.tracers.b))
#for j in 1:grid.Ny
#    b_bcs[j] = buoyancy_relaxation(grid.Nx÷2, j, grid, model.clock, (; b))
#end

fig_1 = Figure(resolution= (600, 800))
ax  = Axis(fig_1[1, 1], title = "Wind Stress Profile", xlabel="τ [Pa]", ylabel="Latitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=-0.1:0.1:0.2, yticks=-70:20:70)
lines = lines!(ax, - τ_bcs .* 1000, φ, linewidth = 3) # (we re-convert the wind stress form kg/m² to Nm)
ylims!(-70, 70)
#display(fig_1)
#save("example/wind_stress_2.png", fig_1)

##########

# Let's plot the initial conditions to make sure they are reasonable
λ = λnodes(cpu_grid.underlying_grid, Center(), Center(), Center())
z = znodes(cpu_grid.underlying_grid, Center(), Center(), Center())

fig = Figure(resolution= (600, 800))
ax  = Axis(fig[1:4, 1], title = "Buoyancy Relaxation Profile",  xlabel="Longitude [∘]", ylabel="Latitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=25, aspect=0.5, xticks=0:20:60, yticks=-70:20:70)
hm  = heatmap!(ax, λ, φ, b[:, :, grid.Nz], colormap = :deep)
cb  = Colorbar(fig[1:4, 2], hm, ticklabelsize = 20, label = L"[m/s^2]", width = 30, labelsize = 20, ticks = 0.0:0.01:0.06)
display(fig)
#save("example/buoyancy_colorbar1.png", fig)



###############
#overlayed plot
#using CairoMakie.AbstractPlotting
#using CairoMakie.AbstractPlotting.MakieLayout
#sing MakieLayout
#using CairoMakie; CairoMakie.activate!()
#using AbstractPlotting

#scene, layout = layoutscene(resolution = (1000, 2000))
#scene = Scene(resolution = (600, 400))
#layout, ax1, ax2 = layout!(scene, 2, 1, widths=[0.5, 0.5])
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
cb = CairoMakie.Colorbar(fig[1, 1], hm, ticklabelsize = 20, label = L"[m/s^2]", labelsize = 20)
# (we re-convert the wind stress form kg/m² to Nm)
f
#save("example/buoyancy_wind_plot_final.png", f)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("buoyancy_wind_plot_final1.png", f, px_per_unit = 3)
#=

using GLMakie, Plots
fig = Figure()
ax2 = Axis(f[1, 1], xticklabelcolor = :red, xtickcolor = :red, xlabel="Longitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=0:20:60, yticks=-70:20:70)
ax1 = Axis(f[1, 1], xticklabelcolor = :blue, xtickcolor = :blue, xlabel="τ [Pa]", ylabel="Latitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=-0.1:0.1:0.2, yticks=-70:20:70)

heatmap!(ax2, λ, φ, b[:, :, grid.Nz], colormap = :thermometer)
lines!(ax1, - τ_bcs .* 1000, φ, linewidth = 3, color = :black) # (we re-convert the wind stress form kg/m² to Nm)

twinx!(ax1)
ax1.xlabel1 = "test"
display(fig)



using CairoMakie
using CairoMakie.AbstractPlotting
using AbstractPlotting
scene = scene(resolution = (600, 400))
layout, ax1, ax2 = layoutscene(scene, 2, 1, split = true)

ax2 = Axis(f[1, 1], xticklabelcolor = :red, xtickcolor = :red, xlabel="Longitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=0:20:60, yticks=-70:20:70)
ax1 = Axis(f[1, 1], xticklabelcolor = :blue, xtickcolor = :blue, xlabel="τ [Pa]", ylabel="Latitude [∘]", xticklabelsize = 20, ylabelsize = 20, xlabelsize=20, yticklabelsize = 20, titlesize=20, aspect=0.5, xticks=-0.1:0.1:0.2, yticks=-70:20:70)

hidexdecorations!(ax2)
hidespines!(ax2)
xaxis_right!(ax2)
linkxaxes!(ax1, ax2)

heatmap!(ax2, λ, φ, b[:, :, grid.Nz], colormap = :thermometer)
lines!(ax1, - τ_bcs .* 1000, φ, linewidth = 3, color = :blue) # (we re-convert the wind stress form kg/m² to Nm)

scene

=#  


