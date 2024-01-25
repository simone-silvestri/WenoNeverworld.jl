using GLMakie, JLD2, Oceananigans, Statistics

# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage4/WenoNeverworldData/weno_thirtytwo_compressed_iteration245338.jld2", "r")
keys(hfile)

grid = NeverWorldGrid(resolution, lat, lon, η)

# Extract data
halo = 0
#z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
#Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
#lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
#lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
η = hfile["η"]["data"]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u = hfile["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = hfile["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])
 

# Buoyancy
blims = (quantile(b[:], 0.2), maximum(b[:]))
Λ = log(blims[1]/blims[2])
contours = blims[2] * exp.(Λ .*  range(0, 1, 11) )
longitudes = [20, 64, 85, 105]

for (i, lon_index) in enumerate(longitudes)
    fig = Figure(resolution=(2000, 1000))
    ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Buoyancy at Longitude = " * string(lon[lon_index]) * "∘", titlesize=50)
    heatmap!(ax, lat, z, b[lon_index, :, :], colormap=:plasma, colorrange=blims)
    contour!(ax, lat, z, b[lon_index, :, :], color=:black, linewidth=3, levels=contours, labels=true,
             labelsize=30, labelfont=:bold, labelcolor=:black)
    display(fig)
    #save("plotting/weno_half_b_" * string(i) * ".png", fig)
end

# Depth integrated TKE
weighted_tke = @. (u^2 + v^2)
depth_integrated_tke = sum(weighted_tke, dims=3)[:, :, 1]

fig = Figure(resolution=(1000, 2000))
ax = Axis(fig[1, 1], xlabel="longitude [∘]", xlabelsize=30, xticklabelsize=30, ylabel="latitude [∘]", ylabelsize=30, title="TKE ", yticklabelsize=30, titlesize=50)
hm = heatmap!(ax, lon, lat, log10.(depth_integrated_tke .+ eps(1000.0)), colorrange=(-1, 3), aspect_ratio=0.8, colormap=:plasma)
cbar1 = Colorbar(fig[1, 2], hm, width=30, ticksize=30)
display(fig)
#save("plotting/weno_half_tke.png", fig)

# One slice of KE
z_index = 69
u_slice = u[:, :, z_index]
v_slice = v[:, :, z_index]
weighted_tke_slice = @. (u_slice^2 + v_slice^2)

fig = Figure(resolution=(5000, 5000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize=20, xticklabelsize=20, ylabel="Latitude [∘]", ylabelsize=20, title="1/32∘", yticklabelsize=20, titlesize=25, aspect=0.5, yticks=-70:10:70)
hm = heatmap!(ax, lon, lat, log10.(weighted_tke_slice .+ eps(1000.0)), colorrange=(-4, 1), colormap=:plasma)
cbar1 = Colorbar(fig[1, 2], hm, width=30, ticksize=30)
display(fig)
save("plotting/tke_slice_thirtytwo.png", fig)

quantile(log10.(weighted_tke_slice .+ eps(1000.0))[:], 0.999)