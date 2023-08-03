using GLMakie
using JLD2, Oceananigans, Statistics      
using Contour                                             
# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage2/WenoNeverworldData/weno_eighth_checkpoint_iteration985500.jld2", "r")
#says eighth but is really half interpolated from 1/4
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
#dz = oceangrid.underlying_grid. [1:end-halo] #oceangrid["dz"]

##
#rossby number/coriolis force math
#lat_arr = Float64.(collect(lat))  #convert lat to an array
f = 2 * 2π / 86400 * sin.(lat .* π / 180)
relu(x) = max(0, x)
Δz = reshape( Δz , (1,1,69))
avg_z = ( (Δz[:, :, 1:end-1] + Δz[:, :, 2:end]) * 0.5 ) 
N² = (b[:, :, 2:end] - b[:, :, 1:end-1]) ./ avg_z
f = reshape(f, (1,280,1))
integrand = avg_z .* sqrt.(relu.(N²) ) ./ (abs.(f) .* π)
#integrand = dz .* sqt.(relu.(dz .* b) ./ (abs.(f) .* π))
integral = sum(integrand, dims=3)[:,:,1]

println(size(integral))
println(size(lat))
println(size(lon))

#reshape to 2D
m, n = size(integral)
ross_2d = integral # reshape(integral, (n, m))
#ross_2d = reshape(integral, (size(lat)[1], size(lon)[1]))


##
fig = Figure(resolution = (1000, 4000))
ax = Axis(fig[1, 1], xlabel="longitude [∘]", xlabelsize = 30, xticklabelsize = 30, ylabel="latitude [∘]", ylabelsize = 30,title="Deformation radius, 1/2∘ " ,  yticklabelsize = 30, titlesize=30)
hm = heatmap!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), colorrange = (3, 5), aspect_ratio = 0.8,colormap = :plasma)

#contour lines
contour_levels = 3
contour_labels = range(4.2, stop=4.8, length=contour_levels)
contour_lines = contour!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), levels=contour_labels, linewidth=3, linecolor=:black, color=:black, labels=true, labelsize=20, labelfont=:bold, labelcolor=:black, labelformatter=(x -> "$(round(10^x/1e3, digits=0))"), constrain_labels=false)


cbar1 = Colorbar(fig[1,2], hm, width = 30, ticksize = 30)
display(fig)
save("plotting/weno_ross_half.png", fig)
##

lat_index = argmin(abs.(-50 .-lat))
lon_index = argmin(abs.(30 .-lon))
ross_2d[lon_index, lat_index]