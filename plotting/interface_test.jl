"""using GLMakie
using JLD2, Oceananigans, Statistics
using ColorSchemes

# Function to create color-coded contour lines for a given resolution
function create_contour_lines!(ax, b, lat, z, contours, color, colormap, label)
    contour!(ax, lat, z, b, levels=contours, linewidth=2, linecolor=color, colormap=colormap, label=label)
end

# Load the data and define buoyancy contours for each resolution
function load_data_and_contours(path::String)
    hfile = jldopen(path)
    oceangrid = hfile["grid"]
    b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
    lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
    z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
    
    blims = (quantile(b[:], 0.2), maximum(b[:]))
    Λ = log(blims[1] / blims[2])
    contours = blims[2] * exp.(Λ .* range(0, 1, 6))
    
    return b, lat, z, contours
end

# Define the file paths for each degree resolution
path_1 = "plotting/weno_eight.jld2" #"/storage2/WenoNeverworldData/weno_one_checkpoint_iteration9993525.jld2" #plotting/weno_fourth.jld2"#fourth data path         #resolution 1 data
path_2 = "/storage2/WenoNeverworldData/weno_one_checkpoint_iteration9993525.jld2" #/storage2/WenoNeverworldData/weno_eighth_checkpoint_iteration985500.jld2"#half data path  #resolution 2 data
degree_resolution_1 = "1/8"
degree_resolution_2 = "1"
lon_index = 10

# Load data and define buoyancy contours for each resolution
b_1, lat_1, z_1, contours_1 = load_data_and_contours(path_1)
b_2, lat_2, z_2, contours_2 = load_data_and_contours(path_2)

# Create the plot
fig = Figure(resolution=(2000, 1000))
#ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours " * string(degree_resolution_1) * "∘ " * string(degree_resolution_2) * "∘", titlesize=50)
ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours at Longitude = " * string(lon_index) * "∘", titlesize=50)


# Plot color-coded contour lines for 1/4-degree resolution with blue colormap
# create_contour_lines!(ax, b_1[lon_index,:, :], lat_1, z_1, contours_1, :blue, :blues, string(degree_resolution_1)* "∘ ")
cl = contour!(ax,lat_1,  z_1, b_1[lon_index,:, :]; levels =  contours_1, colormap =  :blues)
lines!(ax, [1],[1]; color = :blue, label = string(degree_resolution_1)* "∘ ")
# Plot color-coded contour lines for 1/2-degree resolution with red colormap
# create_contour_lines!(ax, b_2[lon_index, :, :], lat_2, z_2, contours_2, :red, :reds, string(degree_resolution_2) * "∘ ")
c2 = contour!(ax,lat_2,  z_2, b_2[lon_index,:, :]; levels =  contours_2, colormap =  :reds)
lines!(ax, [1],[1]; color = :red, label = string(degree_resolution_2)* "∘ ")

axislegend(ax, position=:lb, framecolor=(:grey, 0.5), patchsize=(40, 40), markersize=20, labelsize=20) #colors=Dict(string(degree_resolution_1) * "∘ " => :blue, string(degree_resolution_2) * "∘ " => :red))


display(fig)
save("plotting/interface_3.png" , fig)# * string(degree_resolution_1) * "_" * string(degree_resolution_2) * ".png", fig)
"""

using GLMakie
using JLD2, Oceananigans, Statistics
using ColorSchemes

# Function to create color-coded contour lines for a given resolution
function create_contour_lines!(ax, b, lat, z, contours, color, colormap, label)
    contour!(ax, lat, z, b, levels=contours, linewidth=2, linecolor=color, colormap=colormap, label=label)
end

# Load the data and define buoyancy contours for each resolution
function load_data_and_contours(path::String, lon)
    hfile = jldopen(path)
    oceangrid = hfile["grid"]
    b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
    lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
    z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
    
    blims = (quantile(b[:], 0.2), maximum(b[:]))
    Λ = log(blims[1] / blims[2])
    contours = blims[2] * exp.(Λ .* range(0, 1, 6))
    
    lon_index = argmin(abs.(lon .- lat)) # Find the closest longitude index
    return b, lat, z, contours, lon_index
end

# Define the file paths for each degree resolution
path_1 = "plotting/weno_eight.jld2"
path_2 = "plotting/weno_fourth.jld2"
degree_resolution_1 = "1/8"
degree_resolution_2 = "1/4"
lon = 30

# Load data and define buoyancy contours for each resolution
b_1, lat_1, z_1, contours_1, lon_index_1 = load_data_and_contours(path_1, lon)
b_2, lat_2, z_2, contours_2, lon_index_2 = load_data_and_contours(path_2, lon)

# Create the plot
fig = Figure(resolution=(2000, 1000))
ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours at Longitude = " * string(lon) * "∘", titlesize=50)

# Plot color-coded contour lines for 1/8-degree resolution with blue colormap
create_contour_lines!(ax, b_1[lon_index_1, :, :], lat_1, z_1, contours_1, :blue, :blues, string(degree_resolution_1) * "∘ ")

# Plot color-coded contour lines for 1/4-degree resolution with red colormap
create_contour_lines!(ax, b_2[lon_index_2, :, :], lat_2, z_2, contours_2, :red, :reds, string(degree_resolution_2) * "∘ ")

axislegend(ax, position=:lb, framecolor=(:grey, 0.5), patchsize=(40, 40), markersize=20, labelsize=20)

display(fig)
save("plotting/interface_$(lon)_6.png", fig)
