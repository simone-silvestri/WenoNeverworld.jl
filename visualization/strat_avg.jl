using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute
using CairoMakie
using JLD2, Statistics, Printf

# Define simulation configurations
simulations = [
    ("half_degree",     "/storage4/WenoNeverworldData/half_degree_new/",     "weno_half_check",     "1/2∘", "weno_half__checkpoint_iteration13154400.jld2"),
    ("quarter_degree",  "/storage4/WenoNeverworldData/quarter_degree_new/",      "weno_quarter__check",  "1/4∘", "weno_quarter__checkpoint_iteration70005600.jld2"),
    ("eighth_degree",   "/storage4/WenoNeverworldData/eighth_degree_new/",       "weno_eighth_check",   "1/8∘", "weno_eighth_checkpoint_iteration28667520.jld2"),
    ("sixteenth_degree","/storage4/WenoNeverworldData/sixteenth_degree/",    "weno_sixteen_check","1/16∘", "weno_sixteen_checkpoint_iteration2417760.jld2")
]

for (name, dir, prefix, title, file) in simulations
    @info "Processing simulation: $name"

    variables = ("u", "v", "w", "b")
    neverworld_fields = all_fieldtimeseries(prefix, dir; checkpointer=true, variables, number_files=25)
    grid = neverworld_fields[:u].grid
    times = neverworld_fields[:u].times
    Nt = length(times)
    iterations = Nt-25:Nt

    u = neverworld_fields[:u]
    v = neverworld_fields[:v]
    w = neverworld_fields[:w]
    b = neverworld_fields[:b]

    ū = Diagnostics.time_average(u, iterations)
    v̄ = Diagnostics.time_average(v, iterations)
    w_avg = Diagnostics.time_average(w, iterations)
    b̄ = Diagnostics.time_average(b, iterations)

    # Load the corresponding checkpoint file for lat/z info
    checkpoint_file = joinpath(dir, file)
    @info "Opening checkpoint file: $checkpoint_file"
    hfile = jldopen(checkpoint_file, "r")
    oceangrid = hfile["grid"]
    halo = 7
    z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
    lat = collect(oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo])

    # Compute zonal average of buoyancy
    b_r = mean(interior(b), dims=4)[:,:,:,1] * 100
    b_zonal = dropdims(mean(b_r, dims=1); dims=1)

    # Plot setup
    blims = (quantile(b_r[:], 0.2), maximum(b_r[:]))
    contours_log = [0.5, 0.75, 1, 1.5, 2, 3, 4, 5]

    fig = Figure(resolution=(1000, 2000))
    ax = Axis(fig[1, 1],
              xlabel="Latitude [∘]", ylabel="Depth [m]",
              xlabelsize=20, ylabelsize=20,
              yticks=-4000:1000:0, xticklabelsize=20,
              yticklabelsize=20, xticks=-70:10:70,
              title=title, titlesize=20, aspect=2.0)

    hm = heatmap!(ax, lat, z, b_zonal, colormap=:plasma)
    cont = contour!(ax, lat, z, b_zonal, color=:black, linewidth=3,
             levels=contours_log, labels=false)
    cbar = Colorbar(fig[1, 2], hm, width=25, ticksize=20, ticklabelsize=20,
             label=L"m s^{-2}", labelsize=20, height=Relative(1/5))

    outpath = "figures/zonal_avg_$(name).png"
    save(outpath, fig)
    @info "Saved figure to $outpath"
end
