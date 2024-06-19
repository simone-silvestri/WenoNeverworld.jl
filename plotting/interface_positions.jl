using GLMakie
using JLD2, Oceananigans, Statistics
using ColorSchemes
using LaTeXStrings
using WenoNeverworld
using WenoNeverworld.Diagnostics

# Function to create color-coded contour lines
function create_contour_lines!(ax, b, lat, z, contours, color, colormap, label)
    contour!(ax, lat, z, b, levels=contours, linewidth=50, linecolor=color, colormap=colormap, label=label, labelfont = :bold)
end
halo = 7

# Load data and define buoyancy contours
function load_data_and_contours(path::String)
    hfile = jldopen(path)
    oceangrid = hfile["grid"]
    b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
    lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
    z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
    
    blims = (quantile(b[:], 0.2), maximum(b[:]))
    Λ = log(blims[1] / blims[2])
    contours = blims[2] * exp.(Λ .* range(0, 1, 6))

    #lon_index = argmin(abs.(lon .- lat)) # Find the closest longitude index
    lon_index = argmin(abs.(lon .- oceangrid.underlying_grid.λᶜᵃᵃ))
    #lon_index = argmin(abs.(lon .- oceangrid.underlying_grid.longitude)) # Find the closest longitude index
    return b, lat, z, contours, lon_index
    
end

# Define the file paths
path_1 = "/storage3/WenoNeverworldData/weno_sixteenth_checkpoint_iteration972420.jld2" # #plotting/weno_fourth.jld2"#fourth data path         #resolution 1 data
path_2 =  "/storage2/WenoNeverworldData/weno_one_checkpoint_iteration8106971.jld2" #weno_eighth_checkpoint_iteration985500.jld2" #half data path  #resolution 2 data
degree_resolution_1 = "1/16"
degree_resolution_2 = "1"
lon = 30


# Load data and define buoyancy contours for each resolution
b_1, lat_1, z_1, contours_1, lon_index_1 = load_data_and_contours(path_1)
b_2, lat_2, z_2, contours_2, lon_index_2 = load_data_and_contours(path_2)

#lon_index = argmin(abs.(30 .-lon))

# Create the plot
fig = Figure(resolution=(1000, 2000))
#ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours " * string(degree_resolution_1) * "∘ " * string(degree_resolution_2) * "∘", titlesize=50)
ax = Axis(fig[1, 1], xlabel="Latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="z [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title=L"Interface positions at λ = " * string(lon) * "∘", titlesize=25)


# Plot color-coded contour lines for 1/4-degree resolution with blue colormap
# create_contour_lines!(ax, b_1[lon_index,:, :], lat_1, z_1, contours_1, :blue, :blues, string(degree_resolution_1)* "∘ ")
cl = contour!(ax,lat_1,  z_1, b_1[lon_index_1,:, :]; levels =  contours_1, color = (:darkblue, 0.5), labelfont = :bold, linewidth = 6)
lines!(ax, [1],[1]; color = :darkblue, label = string(degree_resolution_1)* "∘ ")
# Plot color-coded contour lines for 1/2-degree resolution with red colormap
# create_contour_lines!(ax, b_2[lon_index, :, :], lat_2, z_2, contours_2, :red, :reds, string(degree_resolution_2) * "∘ ")
c2 = contour!(ax,lat_2,  z_2, b_2[lon_index_2,:, :]; levels =  contours_2, color = (:red4, 0.5), labelfont = :bold,  linewidth = 6)
lines!(ax, [1],[1]; color = :red4, label = string(degree_resolution_2)* "∘")

axislegend(ax, position=:lb, framecolor=(:grey, 0.5), patchsize=(40, 40), markersize=20, labelsize=20) #colors=Dict(string(degree_resolution_1) * "∘ " => :blue, string(degree_resolution_2) * "∘ " => :red))


display(fig)
save("plotting/interface_$(lon)_7.png" , fig)# * string(degree_resolution_1) * "_" * string(degree_resolution_2) * ".png", fig)