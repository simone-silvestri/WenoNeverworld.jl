using CairoMakie
using JLD2, Oceananigans, Statistics
using LaTeXStrings

# Resolution label => (file path, grid size in x)
resolutions = Dict(
    "1/2°"   => ("/storage4/WenoNeverworldData/half_degree_new/weno_half_checkpoint_iteration42167520.jld2", 280),
    "1/4°"   => ("/storage4/WenoNeverworldData/quarter_degree_new/weno_quarter__checkpoint_iteration70005600.jld2", 560),
    "1/8°"   => ("/storage4/WenoNeverworldData/eighth_degree_new/weno_eighth_checkpoint_iteration28667520.jld2", 1120),
    "1/16°"  => ("/storage4/WenoNeverworldData/sixteenth_degree/weno_sixteen_checkpoint_iteration1156320.jld2", 2240)
)

lines = Dict()  # To store latitude and zonal-mean Rossby radius

lat_cutoff = 5.0

for (label, (filepath, nx)) in resolutions
    @info "Processing $label resolution"
    hfile = jldopen(filepath, "r")

    # Load grid and fields
    oceangrid = hfile["grid"]
    halo = 7
    z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
    Δz = reshape(oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo], 1, 1, :)
    lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
    lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
    b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
    close(hfile)

    # Compute Coriolis frequency
    f = 2 * 2π / 86400 .* sin.(lat .* π / 180)
    f = reshape(f, 1, nx, 1)

    # Compute Brunt–Väisälä frequency and Rossby deformation radius
    avg_z = 0.5 .* (Δz[:, :, 1:end-1] .+ Δz[:, :, 2:end])
    N² = (b[:, :, 2:end] .- b[:, :, 1:end-1]) ./ avg_z
    relu(x) = max(0, x)
    integrand = avg_z .* sqrt.(relu.(N²)) ./ (abs.(f) .* π)
    L_D = sum(integrand, dims=3)[:, :, 1] ./ 1000  # [km]

    # Zonal mean over longitudes
    L_D_zonal_mean = dropdims(mean(L_D, dims=1), dims=1)
    mask = abs.(lat) .> lat_cutoff
    lines[label] = (lat[mask], L_D_zonal_mean[mask])
end

# Plotting all resolutions on the same figure
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1,1], xlabel="Latitude [°]", ylabel = L"L_D \; \mathrm{[km]}",
          xlabelsize=30, ylabelsize=30, titlesize=35,
          xticklabelsize=25, yticklabelsize=25)


#=
for (i, (label, (lat, L_D_lat))) in enumerate(sort(collect(lines); by=first))
    lines!(ax, lat, L_D_lat, label=label, linewidth=4, color=colors[i])
end
=#
# Define consistent styles in correct order
res_order = ["1/2°", "1/4°", "1/8°", "1/16°"]
colors = [:red3, :orange, :green, :purple]

for (i, label) in enumerate(res_order)
    lat, L_D_lat = lines[label]

    # Split into hemispheres
    south_mask = lat .< -lat_cutoff
    north_mask = lat .> lat_cutoff

    if any(south_mask)
        lines!(ax, lat[south_mask], L_D_lat[south_mask],
            linewidth=1.5, color=colors[i])
    end
    if any(north_mask)
        lines!(ax, lat[north_mask], L_D_lat[north_mask],
               label=label, linewidth=1.5, color=colors[i])
    end
end


axislegend(ax; position=:rt, labelsize=25)
fig
save("figures/rossby_radius_vs_latitude.png", fig, px_per_unit=4)


