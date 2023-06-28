using Pkg
using Oceananigans
using JLD2
using CairoMakie
using Printf

filename = "uq_one_degree_salinity_snapshots.jld2" 
file = jldopen(filename)
vars = keys(file["timeseries"]) 
iterations = parse.(Int,keys(file["timeseries/t"]))
last_index = length(iterations)

u_series = FieldTimeSeries(filename, "u"; architecture = CPU())
v_series = FieldTimeSeries(filename, "v"; architecture = CPU())
w_series = FieldTimeSeries(filename, "w"; architecture = CPU())
T_series = FieldTimeSeries(filename, "T"; architecture = CPU())
S_series = FieldTimeSeries(filename, "S"; architecture = CPU())

times = u_series.times

xu, yu, zu = nodes(u_series[1])
xv, yv, zv = nodes(v_series[1])
xT, yT, zT = nodes(T_series[1])

u = u_series[last_index]
v = v_series[last_index]
T = T_series[last_index]

# Compute the derivative of v.

dv = compute!(Field(∂x(v)))
dv_interior = interior(dv, :, :, size(dv, 3)) 

# Compute the integral of temperature in the vertical direction.

int_T_z = compute!(Field(Integral(T, dims = 3)))
int_T_z_interior = interior(int_T_z, :, :, size(int_T_z, 3)) 

# Compute the integral of temperature in all directions.

int_T_xyz = compute!(Field(Integral(T)))

@printf("The heat content over the entire domain is %.6g.", int_T_xyz[1,1,1])

# Plot the surface of dv.
fig = Figure(resolution = (1000,750))

ax_T = Axis(fig[1,1]; xlabel = "x (m)", ylabel = " y (m)", aspect = 1.0, title = "Heat Content")
hm_T = heatmap!(ax_T, xT, yT, int_T_z_interior)
Colorbar(fig[1,2], hm_T)

save("HeatContent.pdf",fig)

nothing # hide