using JLD2, Oceananigans, Statistics, CairoMakie

# Setup
halo = 7
nfiles = 25  # number of most recent files to average

res_paths = Dict("1/16" => "/storage4/WenoNeverworldData/sixteenth_degree/",
                 "1/8"  => "/storage4/WenoNeverworldData/eighth_degree_interp/")

# Get sorted list of latest checkpoint files
function sorted_checkpoint_files(path::String; nfiles::Int)
    all_files = readdir(path; join=true)
    checkpoint_files = filter(f -> occursin("checkpoint", f) && endswith(f, ".jld2"), all_files)

    if isempty(checkpoint_files)
        error("No checkpoint .jld2 files found in $path")
    end

    # Sort by filename (natural order), or use by mod time if needed
    sorted = sort(checkpoint_files)  # or use sort(..., by = f -> stat(f).mtime) for time-based sort
    return sorted[max(end - nfiles + 1, 1):end]
end

# Compute time-averaged depth-integrated KE
function time_averaged_depth_integrated_ke(res::String)
    path = res_paths[res]
    files = sorted_checkpoint_files(path; nfiles=nfiles)

    total_tke = nothing
    lon, lat = nothing, nothing  # Declare here to avoid scope issues

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

        # Staggered v average
        v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])

        ke = @. (u^2 + v^2) * Δz
        depth_integrated = sum(ke, dims=3)[:, :, 1]

        if isnothing(total_tke)
            total_tke = depth_integrated
        else
            total_tke .+= depth_integrated
        end

        close(file)
    end

    return total_tke ./ length(files), lon, lat
end

# Main execution
@info "Computing time-averaged depth-integrated KE for both resolutions..."
ke_16, lon16, lat = time_averaged_depth_integrated_ke("1/16")
ke_8, lon, lat = time_averaged_depth_integrated_ke("1/8")


ke_16_downsampled = ke_16[1:2:end, 1:2:end]
lon_downsampled = lon16[1:2:end]
lat_downsampled = lat[1:2:end]
##### change this depending on the simulation

# Now compute the difference
@info "Computing difference..."
ke_diff = ke_16_downsampled .- ke_8


@info "Plotting difference..."
fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize=40, xticklabelsize=40, title = "1/16° - 1/8°",
          ylabel="Latitude [∘]", ylabelsize=40,
          yticklabelsize=40, titlesize=45, aspect=0.5, yticks=-70:10:70,
          yticksize=15, xticksize=15)

hm = heatmap!(ax, lon, lat, ke_diff, colormap=:balance, colorrange=(-1e3, 1e3))
Colorbar(fig[1,2], hm, width=60, ticksize=40, ticklabelsize=40,
         label=L"m^2/s^2", labelsize=40, height=Relative(3/5))

display(fig)
save("figures/TKE_diff_8_interp.png", fig, px_per_unit=5)
