using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using CairoMakie
using JLD2, Statistics, Printf


# 1) your four sims: key, directory, prefix, title, checkpoint file
sims = [
  #("half_degree",      "/storage4/WenoNeverworldData/half_degree_new/",    "weno_half_check",    "1/2∘",  "weno_half__checkpoint_iteration13154400.jld2"),
  #("quarter_degree",   "/storage4/WenoNeverworldData/quarter_degree_new/","weno_quarter__check", "1/4∘",  "weno_quarter__checkpoint_iteration70005600.jld2"),
  ("eighth_degree",    "/storage4/WenoNeverworldData/eighth_degree_new/", "weno_eighth_check",   "1/8∘",  "weno_eighth_checkpoint_iteration28667520.jld2"),
  ("eighth_degree_interp",    "/storage4/WenoNeverworldData/eighth_degree_interp/", "weno_eighth_check",   "1/8∘ - int",  "weno_eighth_checkpoint_iteration4204801.jld2"),
  #("sixteenth_degree", "/storage4/WenoNeverworldData/sixteenth_degree/",  "weno_sixteen_check",  "1/16∘","weno_sixteen_checkpoint_iteration2417760.jld2")
]

# 2) load each, average last 25, extract b_zonal + grid
b_zonal = Dict{String,Matrix{Float64}}()
# 0) Define halo and load lat/z up front from highest resolution
halo = 7
grid = jldopen("/storage4/WenoNeverworldData/eighth_degree_new/weno_eighth_checkpoint_iteration28667520.jld2", "r") do h
    h["grid"]
end
lat = collect(grid.underlying_grid.φᵃᶜᵃ[1:end-halo])
z   = grid.underlying_grid.zᵃᵃᶜ[1:end-halo]


for (key, dir, prefix, title, file) in sims
  @info "Processing $title"

    fts = all_fieldtimeseries(prefix, dir;
            checkpointer=true,
            variables=("u","v","w","b"),
            number_files=1)
    if fts[:b] === nothing
        @warn "No field timeseries found for $key. Skipping..."
        continue
    end

  Nt = length(fts[:u].times)
  its = Nt-0:Nt
  b̄  = Diagnostics.time_average(fts[:b], its)

  ck = jldopen(joinpath(dir,file),"r") do h
    h["grid"]
  end

  b_r          = mean(interior(b̄), dims=4)[:,:,:,1] * 100
  b_zonal[key] = dropdims(mean(b_r, dims=1); dims=1)
end

# 3) Sanity-check
@info "Sizes:"
for key in keys(b_zonal)
  @info "  $key => $(size(b_zonal[key]))"
end
#=

# 4) Ensure downsampling matches corresponding field sizes
b16_8 = b16[1:2:end, :]
lat_8 = lat[1:2:end]  # downsample latitudes by 2
lat_8 = lat_8[1:size(b16_8, 1)]  # Ensure latitudes match the size of b16_8

b16_4 = b16[1:4:end, :]
lat_4 = lat[1:4:end]  # downsample latitudes by 4
lat_4 = lat_4[1:size(b16_4, 1)]  # Ensure latitudes match the size of b16_4

b16_2 = b16[1:8:end, :]
lat_2 = lat[1:8:end]  # downsample latitudes by 8
lat_2 = lat_2[1:size(b16_2, 1)]  # Ensure latitudes match the size of b16_2
=#



# 5) Compute diffs on their matching lat grids
diffs = Dict(
   #"1/16−1/8" => (b16_8 .- b_zonal["eighth_degree"],   lat_8),
   #1/16−1/4" => (b16_4 .- b_zonal["quarter_degree"],  lat_4),
   #"1/16−1/2" => (b16_2 .- b_zonal["half_degree"],     lat_2),
   "1/8-1/8 - int" => (b_zonal["eighth_degree"] - b_zonal["eighth_degree_interp"], lat)
   )


levels = [-1, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]

xticks = (-70:10:70, string.(-70:10:70))

# 6) Plot & save each
for (tag, (diff, latsub)) in diffs
  fig = Figure(resolution=(1000, 2000))
  ax  = Axis(fig[1,1],
             title="1/8° - 1/8° (interpolated)",
             xlabel="Latitude (∘)", ylabel="Depth (m)", xlabelsize = 20, ylabelsize =20,
             xticks=xticks, yticks=-4000:1000:0, xticklabelsize=20, yticklabelsize=20, 
             titlesize=25, 
             aspect=2)


  hm = heatmap!(ax, latsub, z, diff, colormap=:plasma)
  contour!(ax, latsub, z, diff, levels=levels, color=:black, linewidth=2)
  text!(ax, "70", position = (70, -3000), align = (:center, :top), fontsize = 20)
  Colorbar(fig[1,2], hm,
           label="m s⁻²", labelsize=20, ticklabelsize=16,
           width=25, height=Relative(1/5))

safe_tag = replace(tag, "_" => "-")  # Make sure no special characters are in the file name
out = "figures/zonal_avg_diff_8interp.png"
  save(out, fig)
  @info "Saved $out"
end



