using JLD2, Oceananigans, Statistics, CairoMakie

# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage4/WenoNeverworldData/sixteenth_degree_new/weno_sixteen_checkpoint_iteration1156320.jld2", "r")
keys(hfile)
#initialized from 1/4
## grab grid and fields
oceangrid = hfile["grid"]
## 
halo = 7
z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u = hfile["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = hfile["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])

#buoyancy
blims = (quantile(b[:], 0.2), maximum(b[:]))
Λ = log(blims[1]/blims[2])
contours = blims[2] * exp.(Λ .*  range(0, 1, 11) )
longitudes = [160, 512, 680, 840]
for (i,lon_index) in enumerate(longitudes) 
    fig = Figure(resolution = (2000, 1000))
    ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize = 30, xticks = -60:10:60, xticklabelsize = 30, ylabel="depth [m]", ylabelsize = 30, yticks = -4000:1000:0, yticklabelsize = 30, title="Buoyancy at Longitude = " * string(lon[lon_index]) * "∘", titlesize=50)
    heatmap!(ax, lat, z, b[lon_index, :, :], colormap= :plasma,  colorrange = blims)
    contour!(ax, lat, z, b[lon_index, :, :], color=:black, linewidth=3, levels=contours, labels = true,
    labelsize = 30, labelfont = :bold, labelcolor = :black)
    #display(fig)
    #save("plotting/weno_half_b_" * string(i) * ".png", fig)
end


#########
# Depth integrated TKE 

m, n, ℓ = size(u)
Δz = reshape(Δz,  (1,1,ℓ))
weighted_tke = @. (u^2 + v^2) * Δz
depth_integrated_tke = sum(weighted_tke, dims=3)[:,:,1]

fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/16∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:10:70, yticksize = 15, xticksize = 15)
hm = heatmap!(ax, lon, lat, log10.(depth_integrated_tke .+ eps(1000.0)), colorrange = (-1, 3), colormap = :plasma)
cbar1 = Colorbar(fig[1,2], hm, width = 60, ticksize = 40, ticklabelsize = 40, label=L"m^{2} s^{-2}", labelsize = 40, height = Relative(2/3))
display(fig)
CairoMakie.activate!()
CairoMakie.save("figures/TKE_sixteen_cb.png", fig, px_per_unit = 5)
#########


#mean(weighted_tke_slice)
#=
##
#One slice of KE
m, n, ℓ = size(u)
Δz = reshape(Δz,  (1,1,ℓ))

#weighted_tke = (u^2 + v^2) * Δz
z_index = 34
u_slice = u[:, :, z_index]
v_slice = v[:, :, z_index]
weighted_tke_slice = (u_slice .^ 2 + v_slice .^ 2)

##
fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], aspect = 0.5, xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/16∘", yticklabelsize = 40, titlesize=45, yticks=-70:10:70, yticksize = 15, xticksize = 15)
hm = heatmap!(ax, lon, lat, log10.(weighted_tke_slice .+ eps(1000.0)), colorrange = (-4, 1), colormap = :plasma)
cbar1 = Colorbar(fig[1,2], hm, width = 60, ticksize = 40, ticklabelsize = 40, label=L"m^{2} s^{-2}", labelsize = 40, height = Relative(2/3))
display(fig)

using CairoMakie
CairoMakie.activate!()
CairoMakie.save("figures/KE_sixteen_cb.png", fig, px_per_unit = 5)
##
quantile(log10.(weighted_tke_slice .+ eps(1000.0))[:], 0.999)
=#