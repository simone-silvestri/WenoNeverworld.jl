using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute
using JLD2, Oceananigans, Statistics
using GLMakie
#sing CairoMakie

dir = "../"
# stride = 3

#Load the data
@info "Loading data..."
local_path = pwd()
#add the correct file to use as a grid
hfile_1 = jldopen("/storage4/WenoNeverworldData/weno_thirtytwo_compressed_iteration_new1269929.jld2", "r")
keys(hfile_1)
#oceangrid_1 = hfile_1["grid"]

resolution = 1/32
oceangrid_1 = NeverworldGrid(resolution)
halo = 7
z = oceangrid_1.underlying_grid.zᵃᵃᶜ[1:end-halo]
lat = collect(oceangrid_1.underlying_grid.φᵃᶜᵃ[1:end-halo])

using HDF5
hfile = h5open("plotting/b_avg_data.h5", "r")
b_avgs = read(hfile["b_avg"])
print(keys(hfile))
#b_avgs = hfile["b_avg"]

#b_r = mean(b_avgs, dims = 4)[:,:,:,1] * 100
b_r = b_avgs
blims = (quantile(b_r[:], 0.2), maximum(b_r[:]))


Λ = log(blims[1]/blims[2])

contours_log = [0.5, 0.75, 1, 1.5, 2, 3, 4, 5]

lon_index = round(Int, size(b_avgs)[1]/2)
fig = Figure(resolution=(1000,2000))
ax = Axis(fig[1, 1], xlabel="Latitude [∘]", xlabelsize=30, yticks=-5000:1000:0, xticklabelsize=30, ylabel="Depth [m]", ylabelsize=30, xticks=-90:20:1120, yticklabelsize=30, title="1/32∘", titlesize=50, aspect=2.0)

print(size(lat))
print(size(z))
#print(size(b_r[lon_index, :, :]))

hm = GLMakie.heatmap!(ax, lat, z, b_r[lon_index, :, :], colormap= :plasma, levels =9)  #colorrange = (7.5e-5, 3.6e-5)
#hm = CairoMakie.heatmap!(ax, lat[1:2:end], z[1:2:end], b_r[lon_index, 1:2:end, 1:2:end], colormap=:plasma, levels=9)

GLMakie.contour!(ax, lat, z,  b_r[lon_index, :, :], color=:black, linewidth=3, levels=contours_log, labels = true,labelsize = 30, labelfont = :bold, labelcolor = :black)

GLMakie.save("plotting/ta_strat_thirtytwo_final_3.png", fig)
close(hfile)