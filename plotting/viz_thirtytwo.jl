
using WenoNeverworld
using GLMakie, JLD2, Oceananigans, Statistics

# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage4/WenoNeverworldData/weno_thirtytwo_compressed_iteration_new75647.jld2", "r")
keys(hfile)

resolution = 1/32

oceangrid = NeverworldGrid(resolution)

halo = 0
z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u = hfile["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = hfile["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])


#mean(weighted_tke_slice)

##
#One slice of KE
m, n, ℓ = size(u)


#weighted_tke = (u^2 + v^2) * Δz
# z_index = level => 69 = surface slice
z_index = 69
u_slice = u[:, :, z_index]
v_slice = v[:, :, z_index]
weighted_tke_slice = (u_slice .^ 2 + v_slice .^ 2)
#weighted_tke_slice = (u_slice .^ 2 + v_slice .^ 2) .* Δz[:, :, z_index]

##
fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/32∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:10:70, yticksize = 15, xticksize = 15)
hm = heatmap!(ax, lon, lat, log10.(weighted_tke_slice .+ eps(1000.0)), colorrange = (-4, 1), colormap = :plasma)
#cbar1 = Colorbar(fig[1,2], hm, width = 30, ticksize = 30, label=L"[m$^2$/s$^2$]", labelsize = 20, ticklabelsize = 20)
display(fig)
#save("plotting/tke_slice_thirtytwo.png", fig)
##
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("KE_thirtytwo_new.png", fig, px_per_unit = 5)

quantile(log10.(weighted_tke_slice .+ eps(1000.0))[:], 0.999)