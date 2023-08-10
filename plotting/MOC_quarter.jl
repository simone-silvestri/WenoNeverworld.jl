using GLMakie
using JLD2, Oceananigans, Statistics      
                                            
# Load the data
@info "Loading data..."
local_path = pwd()
hfile = jldopen("/storage2/WenoNeverworldData/weno_one_checkpoint_iteration3764471.jld2", "r")
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
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u = hfile["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = hfile["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])
degree_resolution = "1"
lon = 30
level = 30
# Create the plot
fig = Figure(resolution=(2000, 1000))
#ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours " * string(degree_resolution_1) * "∘ " * string(degree_resolution_2) * "∘", titlesize=50)
ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="Interface Contours at Longitude = " * string(lon) * "∘" * string(degree_resolution) * "∘" , titlesize=50)

# MOC math
function calculate_eulerian_MOC(v::AbstractArray)
    v̄ = mean(v, dims=1)
    ψ = similar(v)

    for k in 2:size(v, 3)
        dz = Δz[k - 1]
        for j in 1:size(v, 2)
            for i in 1:size(v, 1)
                ψ[i, j, k] = ψ[i, j, k - 1] + v̄[1, j, k - 1] * dz
            end
        end
    end

    return ψ
end
##
moc_field = calculate_eulerian_MOC(v)  

GLMakie.contour!(ax, lat, z, moc_field[level,:,:])


display(fig)
save("plotting/moc_one_$(lon)_.png", fig)

