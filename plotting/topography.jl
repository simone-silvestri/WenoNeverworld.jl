using GLMakie, JLD2, Oceananigans, Statistics
using WenoNeverworld
# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage3/WenoNeverworldData/weno_sixteenth_checkpoint_iteration972420.jld2", "r")
keys(hfile)
#initialized from 1/4
## grab grid and fields
oceangrid = hfile["grid"]
## 
halo = 7
#z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]

bathymetry_params = WenoNeverworld.NeverWorldBathymetryParameters(; scotia_arc = WenoNeverworld.ScotiaArcParameters(; depth = 3000)) 
grid = NeverworldGrid(1/16; bathymetry_params)
z = interior(grid.immersed_boundary.bottom_height, :, :, 1)

print(size(z), size(lat), size(lon))

##
fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 20, xticklabelsize = 20, ylabel="Latitude [∘]", ylabelsize = 20,title="Topography",  yticklabelsize = 20, titlesize=25, aspect=0.5, yticks=-70:20:70, xticks=0:20:60, background_color = :transparent)
hm = heatmap!(ax, lon, lat, z, colorrange = (-4000, 0), colormap = Reverse(:deep))
cbar1 = Colorbar(fig[1,2], hm, width = 30, ticksize = 10, label="Depth [m]", labelsize = 20, ticklabelsize = 20)
display(fig)
save("topo_sixteen_2.png", fig)
##