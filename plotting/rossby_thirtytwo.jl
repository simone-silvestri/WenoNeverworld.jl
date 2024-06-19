using GLMakie
using JLD2, Oceananigans, Statistics      
using WenoNeverworld
#using Contour                                             
# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage4/WenoNeverworldData/weno_thirtytwo_compressed_iteration_new102719.jld2", "r")

keys(hfile)
resolution = 1/32

oceangrid = NeverworldGrid(resolution)
## grab grid and fields
#oceangrid = hfile["grid"]
## 
halo = 7
z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]

halo = 0
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
#dz = oceangrid.underlying_grid. [1:end-halo] #oceangrid["dz"]

##
#rossby number/coriolis force math
#lat_arr = Float64.(collect(lat))  #convert lat to an array
f = 2 * 2π / 86400 * sin.(lat .* π / 180)
relu(x) = max(0, x)
Δz = reshape( Δz , (1,1,69))
print(size(Δz), size(z), size(b))

avg_z = ( (Δz[:, :, 1:end-1] + Δz[:, :, 2:end]) * 0.5 ) 
N² = (b[:, :, 2:end] - b[:, :, 1:end-1]) ./ avg_z
f = reshape(f, (1,4480,1))
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
#plot log def radius with log color bar
fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 20, xticklabelsize = 20, ylabel="Latitude [∘]", ylabelsize = 20,title="1/32∘" ,  yticklabelsize = 20, titlesize=25, aspect=0.5, xticks=0:20:60, yticks=-70:20:70)
hm = heatmap!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), colorrange = (3, 5), aspect_ratio = 0.5,colormap = :plasma)

#contour lines
contour_levels = 3
contour_labels = range(4.2, stop=4.8, length=contour_levels)
contour_lines = GLMakie.contour!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), levels=contour_labels, linewidth=2, linecolor=:black, color=:black, labels=false, labelsize=20, labelfont=:bold, labelcolor=:black, labelformatter=(x -> "$(round(10^x/1e3, digits=0))"), constrain_labels=true)

cbar1 = Colorbar(fig[1,2], hm, width = 30, ticklabelsize = 20, label = "log [deformation radius]", labelsize = 20)
display(fig)
save("ross_thirtytwo_log.png", fig)




##plot def radius with log contours
fig = Figure(resolution = (500, 500))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 20, xticklabelsize = 20, ylabel="Latitude [∘]", ylabelsize = 20,title="1/32∘" ,  yticklabelsize = 20, titlesize=25, aspect=0.5, xticks=0:20:60, yticks=-70:20:70)
#hm = heatmap!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), colorrange = (3, 5), aspect_ratio = 0.5,colormap = :plasma)
hm = heatmap!(ax, collect(lon), collect(lat), ross_2d, aspect_ratio = 0.5,colormap = :plasma, colorrange = (10^3, 10^5))

#contour lines that are log10 of the def radius
contour_levels = 3
contour_labels = range(4.2, stop=4.8, length=contour_levels)
contour_lines = GLMakie.contour!(ax, collect(lon), collect(lat), log10.(ross_2d .+ eps(1000.0)), levels=contour_labels, linewidth=2, linecolor=:black, color=:black, labels=false, labelsize=20, labelfont=:bold, labelcolor=:black, labelformatter=(x -> "$(round(10^x/1e3, digits=0))"), constrain_labels=true)

cbar1 = Colorbar(fig[1,2], hm, width = 30, ticklabelsize = 20, label = "Deformation radius [m]", labelsize = 20, scale = log10)
display(fig)
save("ross_thirtytwo.png", fig)
##

lat_index = argmin(abs.(-50 .-lat))
lon_index = argmin(abs.(30 .-lon))
ross_2d[lon_index, lat_index]