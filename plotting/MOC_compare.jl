using GLMakie
using JLD2, Oceananigans, Statistics      

# Load the data for the half-degree resolution
@info "Loading data..."
local_path = pwd()


lon = 30
level = 5


hfile_1 = jldopen("/storage2/WenoNeverworldData/weno_one_checkpoint_iteration5212104.jld2", "r")
oceangrid_1 = hfile_1["grid"]
halo = 7
z_1 = oceangrid_1.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz_1 = oceangrid_1.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat_1 = collect(oceangrid_1.underlying_grid.φᵃᶜᵃ[1:end-halo])
#lon_1 = collect(oceangrid_1.underlying_grid.λᶜᵃᵃ[1:end-halo])
lon_index_1 = argmin(abs.(lon .- oceangrid.underlying_grid.λᶜᵃᵃ))
η_1 = hfile_1["η"]["data"][halo:end-halo, halo:end-halo]
b_1 = hfile_1["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u_1 = hfile_1["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v_1 = hfile_1["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v_1 = 0.5 * (v_1[:, 1:end-1, :] + v_1[:, 2:end, :])
degree_resolution_1 = "1"


hfile_2 = jldopen("/storage2/WenoNeverworldData/weno_half_checkpoint_iteration985500.jld2", "r")   #weno_fourth_checkpoint_iteration689186.jld2", "r") 
oceangrid_2 = hfile_2["grid"]
halo = 7
z_2 = oceangrid_2.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz_2 = oceangrid_2.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat_2 = collect(oceangrid_2.underlying_grid.φᵃᶜᵃ[1:end-halo])
#lon_2 = collect(oceangrid_2.underlying_grid.λᶜᵃᵃ[1:end-halo])
lon_index_2 = argmin(abs.(lon .- oceangrid.underlying_grid.λᶜᵃᵃ))
η_2 = hfile_2["η"]["data"][halo:end-halo, halo:end-halo]
b_2 = hfile_2["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u_2 = hfile_2["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v_2 = hfile_2["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v_2 = 0.5 * (v_2[:, 1:end-1, :] + v_2[:, 2:end, :])
degree_resolution_2 = "1/2"


# Create the plot
fig = Figure(resolution=(2000, 1000))
ax = Axis(fig[1, 1], xlabel="latitude [∘]", xlabelsize=30, xticks=-60:10:60, xticklabelsize=30, ylabel="depth [m]", ylabelsize=30, yticks=-4000:1000:0, yticklabelsize=30, title="MOC Contour Comparison at Longitude = $lon∘", titlesize=50)

# Calculate MOC for degree 1
function calculate_eulerian_MOC_1(v_1::AbstractArray, Δz_1::Vector{Float64})
    v̄_1 = mean(v_1, dims=1)
    ψ_1 = similar(v_1)

    for k in 2:size(v_1, 3)
        dz_1 = Δz_1[k - 1]
        for j in 1:size(v_1, 2)
            for i in 1:size(v_1, 1)
                ψ_1[i, j, k] = ψ_1[i, j, k - 1] + v̄_1[1, j, k - 1] * dz_1
            end
        end
    end

    return ψ_1
end

moc_field_1 = calculate_eulerian_MOC_1(v_1, Δz_1)
##



##
# Calculate MOC for degree 2
function calculate_eulerian_MOC_2(v_2::AbstractArray, Δz_2::Vector{Float64})
    v̄_2 = mean(v_2, dims=1)
    ψ_2 = similar(v_2)

    for k in 2:size(v_2, 3)
        dz_2 = Δz_2[k - 1]
        for j in 1:size(v_2, 2)
            for i in 1:size(v_2, 1)
                ψ_2[i, j, k] = ψ_2[i, j, k - 1] + v̄_2[1, j, k - 1] * dz_2
            end
        end
    end

    return ψ_2
end

##
moc_field_2 = calculate_eulerian_MOC_2(v_2, Δz_2)


GLMakie.contour!(ax, lat_1, z_1, moc_field_1[lon_index_1,:, :],levels = level, color=:blue)
GLMakie.contour!(ax, lat_2, z_2, moc_field_2[lon_index_2,:, :], levels = level, color=:red)
lines!(ax, [1],[1]; color = :blue, label = string(degree_resolution_1)* "∘ ")
lines!(ax, [1],[1]; color = :red, label = string(degree_resolution_2)* "∘ ")

axislegend(ax, position=:lb, framecolor=(:grey, 0.5), patchsize=(40, 40), markersize=20, labelsize=20) #colors=Dict(string(degree_resolution_1) * "∘ " => :blue, string(degree_resolution_2) * "∘ " => :red))


display(fig)
save("plotting/moc_compare_$(lon)_1.png", fig)


