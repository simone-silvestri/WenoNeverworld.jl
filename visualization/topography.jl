using CairoMakie, JLD2, Oceananigans, Statistics
using WenoNeverworld
# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage4/WenoNeverworldData/quarter_degree_new/weno_quarter__checkpoint_iteration70005600.jld2", "r")
keys(hfile)
#initialized from 1/4
## grab grid and fields
oceangrid = hfile["grid"]
## 
halo = 7
#z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]

bathymetry_params = WenoNeverworld.NeverWorldBathymetryParameters() 
grid = NeverworldGrid(1/4; bathymetry_params)
z = interior(grid.immersed_boundary.bottom_height, :, :, 1)

print(size(z), size(lat), size(lon))

##
fig = Figure(resolution = (600, 800))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 20, xticklabelsize = 20, ylabel="Latitude [∘]", ylabelsize = 20,  yticklabelsize = 20, titlesize=25, aspect=0.5, yticks=-70:20:70, xticks=0:20:60, backgroundcolor = :transparent)
hm = heatmap!(ax, lon, lat, z, colorrange = (-3000, 0), colormap = Reverse(:deep))
cbar1 = Colorbar(fig[1,2], hm, width = 30, ticksize = 10, label="Depth [m]", labelsize = 20, ticklabelsize = 20)
display(fig)

CairoMakie.activate!()
CairoMakie.save("topography.png", fig, px_per_unit = 3)
