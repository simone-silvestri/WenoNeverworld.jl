using JLD2, Oceananigans, Statistics, CairoMakie

# Setup
halo = 7
nfiles = 25  # number of most recent files to average

res_paths = Dict(
    "1/16" => "/storage4/WenoNeverworldData/sixteenth_degree/",
    "1/8"  => "/storage4/WenoNeverworldData/eighth_degree_new/",
    "1/4"  => "/storage4/WenoNeverworldData/quarter_degree_new/",
    "1/2"  => "/storage4/WenoNeverworldData/half_degree_new/"
)

# Helper: Sorted checkpoint files
function sorted_checkpoint_files(path::String; nfiles::Int)
    all_files = readdir(path; join=true)
    checkpoint_files = filter(f -> occursin("checkpoint", f) && endswith(f, ".jld2"), all_files)
    if isempty(checkpoint_files)
        error("No checkpoint .jld2 files found in $path")
    end
    sorted = sort(checkpoint_files)
    return sorted[max(end - nfiles + 1, 1):end]
end

# Compute time-averaged depth-integrated KE
function time_averaged_depth_integrated_ke(res::String)
    path = res_paths[res]
    files = sorted_checkpoint_files(path; nfiles=nfiles)

    total_tke = nothing
    lon, lat = nothing, nothing

    for filepath in files
        file = jldopen(filepath, "r")
        oceangrid = file["grid"]
        u = file["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
        v = file["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
        Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
        if isnothing(lon) || isnothing(lat)
            lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
            lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
        end
        Δz = reshape(Δz, 1, 1, :)

        v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])
        ke = @. (u^2 + v^2) * Δz
        depth_integrated = sum(ke, dims=3)[:, :, 1]

        total_tke = isnothing(total_tke) ? depth_integrated : total_tke .+ depth_integrated
        close(file)
    end

    return total_tke ./ length(files), lon, lat
end

# Plotting function
function plot_ke_map(tke, lon, lat, res_label)
    fig = Figure(resolution = (1000, 2000))
    ax = Axis(fig[1, 1],
        xlabel="Longitude [∘]", ylabel="Latitude [∘]",
        xlabelsize=40, ylabelsize=40,
        xticklabelsize=40, yticklabelsize=40,
        titlesize=45, title= res_label * "°",
        aspect=0.5, yticks=-70:10:70, yticksize=15, xticksize=15)

    hm = heatmap!(ax, lon, lat, tke, colormap=:plasma)
    Colorbar(fig[1, 2], hm, width=60, ticksize=40, ticklabelsize=40,
             label=L"m^2/s^2", labelsize=40, height=Relative(3/5))

    display(fig)
    save("figures/TKE_avg_$(replace(res_label, '/' => '_')).png", fig, px_per_unit=5)
end

# Loop through resolutions and plot
for res in keys(res_paths)
    @info "Processing $res resolution..."
    ke, lon, lat = time_averaged_depth_integrated_ke(res)
    plot_ke_map(ke, lon, lat, res)
end
