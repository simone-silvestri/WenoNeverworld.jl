using GLMakie
using JLD2, Oceananigans, Statistics      
using Contour                                             
# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("plotting/weno_fourth.jld2")
#hfile = jldopen("/storage2/WenoNeverworldData/weno_four_checkpoint_iteration995943.jld2", "r")

keys(hfile)
#initialized from 1/4
## grab grid and fields
oceangrid = hfile["grid"]
## 
halo = 7
z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat = collect(oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo])
lon = collect(oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo])
#lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
#lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
#dz = oceangrid.underlying_grid. [1:end-halo] #oceangrid["dz"]


##
#rossby number/coriolis force math
#lat_arr = Float64.(collect(lat))  #convert lat to an array
f = 2 * 2π / 86400 * sin.(lat .* π / 180)
relu(x) = max(0, x)
Δz = reshape( Δz , (1,1,69))
avg_z = ( (Δz[:, :, 1:end-1] + Δz[:, :, 2:end]) * 0.5 ) 
N² = (b[:, :, 2:end] - b[:, :, 1:end-1]) ./ avg_z
f = reshape(f, (1,560,1))
integrand = avg_z .* sqrt.(relu.(N²) ) ./ (abs.(f) .* π)
#integrand = dz .* sqt.(relu.(dz .* b) ./ (abs.(f) .* π))
integral = sum(integrand, dims=3)[:,:,1]

println(size(integral))


#reshape to 2D
m, n = size(integral)
ross_2d = integral # reshape(integral, (n, m))
#ross_2d = reshape(integral, (size(lat)[1], size(lon)[1]))


##
fig = Figure(resolution=(1000, 1000))  # Reduce the resolution here
ax = Axis(fig[1, 1], xlabel="longitude [∘]", xlabelsize=30, xticklabelsize=30, ylabel="latitude [∘]", ylabelsize=30, title="Deformation radius, 1/4∘", yticklabelsize=30, titlesize=30)
hm = heatmap!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), colorrange=(3, 5), aspect_ratio=0.8, colormap=:plasma)

# Contour lines
contour_levels = 5
contour_labels = range(3.6, stop=4.8, length=contour_levels)
contour_lines = contour!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), levels=contour_labels, linewidth=3, linecolor=:black, color=:black)

# Manually position the contour labels to avoid overlap
for (i, level) in enumerate(contour_labels)
    idx = argmin(abs.(log10.(ross_2d .+ eps(1000.0)) .- level))
    x, y = lon[idx], lat[idx]
    text!(ax, x, y, "10^$(round(level, digits=1))", halign=:center, valign=:center, fontsize=20, color=:black, bold=true)
end

cbar1 = Colorbar(fig[1, 2], hm, width=30, ticksize=30)
display(fig)