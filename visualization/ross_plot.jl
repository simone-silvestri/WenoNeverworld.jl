using CairoMakie
using JLD2, Oceananigans, Statistics, FFTW

# Manual mapping of label → (file path, nx size)
resolutions = Dict(
    "1/2°"   => ("/storage4/WenoNeverworldData/half_degree_new/weno_half_checkpoint_iteration42167520.jld2", 280),
    "1/4°"   => ("/storage4/WenoNeverworldData/quarter_degree_new/weno_quarter__checkpoint_iteration70005600.jld2", 560),
    "1/8°"   => ("/storage4/WenoNeverworldData/eighth_degree_new/weno_eighth_checkpoint_iteration28667520.jld2", 1120),
    "1/16°"  => ("/storage4/WenoNeverworldData/sixteenth_degree_new/weno_sixteen_checkpoint_iteration1156320.jld2", 2240)
)

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

    # Compute Coriolis frequency and reshape
    f = 2 * 2π / 86400 .* sin.(lat .* π / 180)
    f = reshape(f, 1, nx, 1)

    # Compute Rossby deformation radius
    avg_z = 0.5 .* (Δz[:, :, 1:end-1] .+ Δz[:, :, 2:end])
    N² = (b[:, :, 2:end] .- b[:, :, 1:end-1]) ./ avg_z
    relu(x) = max(0, x)
    integrand = avg_z .* sqrt.(relu.(N²)) ./ (abs.(f) .* π)
    integral_km = sum(integrand, dims=3)[:, :, 1] ./ 1000  # [km]

    # Plot
    fig = Figure(resolution = (1000, 2000))
    ax = Axis(fig[1, 1],
        xlabel="Longitude [∘]", ylabel="Latitude [∘]",
        xlabelsize=40, ylabelsize=40,
        xticklabelsize=40, yticklabelsize=40,
        titlesize=45, title=label,
        xticks=0:20:60, yticks=-70:20:70,
        backgroundcolor=:transparent, aspect=0.5)

    hm = heatmap!(ax, collect(lon), collect(lat), integral_km,
        colormap = :plasma, colorrange = (10, 100))

    contour_levels = [20, 40, 60]
    contour!(ax, collect(lon), collect(lat), integral_km;
                 levels=contour_levels, linewidth=5, color=:black)
        
    cbar_ticks_km = [20, 40, 60, 80, 100]  # You can change these if needed
    cbar_labels = string.(cbar_ticks_km)
    cbar1 = Colorbar(fig[1,2], hm; width=30, ticklabelsize=40, label=L"L_D \, [km]", labelsize=40,
                         ticks=(cbar_ticks_km, cbar_labels), height=Relative(3/4))
        

        save("rossby_deformation_$(replace(label, "/" => "over")).png", fig)
end
