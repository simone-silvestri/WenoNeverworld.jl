using JLD2 
# 

data_directory = "/storage4/WenoNeverworldData/quarter_SO_to_NO/"
data_directory = "/home/sandre/Repositories/WenoNeverworld.jl/simulations_quarter/"
data_directory = "/home/sandre/Repositories/WenoNeverworld.jl/"
filename = "wfs.jld2" # memo to self, grab snapshots instead of checkpoint
filename = "weno_eighth_checkpoint_iteration2353875.jld2"
jlfile = jldopen(data_directory * filename)
#=
tkeys = keys(jlfile["timeseries"]["b"])
b = jlfile["timeseries"]["b"][tkeys[end]]
u = jlfile["timeseries"]["u"][tkeys[end]]
v = jlfile["timeseries"]["v"][tkeys[end]]
=#
b = jlfile["b"]["data"]
u = jlfile["u"]["data"]
v = jlfile["v"]["data"]
halo = 7
#=
zs = jlfile["grid"]["underlying_grid"]["zᵃᵃᶜ"][1+halo:end-halo]
λs = jlfile["grid"]["underlying_grid"]["λᶜᵃᵃ"][1+halo:end-halo]
φs = jlfile["grid"]["underlying_grid"]["φᵃᶜᵃ"][1+halo:end-halo]
=#

using GLMakie 
fig = Figure() 
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "z")
slider = Slider(fig[2, 1], range = 1:size(b)[1],  startvalue = 180)
lon = slider.value
field =  @lift(b[$lon, :, :])
# heatmap!(ax, φs, zs, field, colormap = :thermal)
# contour!(ax, φs, zs, field, color = :black, levels = 20)
heatmap!(ax, field, colormap = :thermal)
contour!(ax, field, color = :black, levels = 20)
fig