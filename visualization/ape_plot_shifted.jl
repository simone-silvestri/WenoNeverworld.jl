using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2

# --- List of simulations
simulations = [
   ("weno_half_ch", "", "1/2∘", "weno_half_ch"),
   ("weno_quarter__ch", "", "1/4∘", "weno_quarter__ch"),
   ("weno_eighth_checkpoint_iteration", "", "1/8∘", "weno_eighth"),
   ("weno_eighth_checkpoint_iteration", "", "1/8∘ - int", "eighth_interp"),
   ("weno_sixteen_checkpoint_iteration", "", "1/16∘ - int", "sixteenth_interp"),
]

# --- First, load last time from quarter-degree simulation
quarter_filename = "weno_quarter__ch_ape_with_time.jld2"
quarter_data = load(quarter_filename)
last_time_quarter = last(quarter_data["time_years"])

println("Last 1/4° time (years): ", last_time_quarter)

# --- Now, load last time from eighth-degree simulation
eighth_filename = "eighth_interp_ape_with_time.jld2"
eighth_data = load(eighth_filename)
last_time_eighth = last(eighth_data["time_years"])

println("Last 1/8° time (years): ", last_time_eighth)


colors = [:red3, :darkorange, :green, :navy, :purple]


# --- Now plot everything
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
          xlabel = L"Time \, [yrs]",
          ylabel = L"APE \; [m^5\,s^{-2}]",
          xlabelsize = 25,
          ylabelsize = 25,
          xticklabelsize = 20,
          yticklabelsize = 20,  
          limits = (nothing, (-7.750e22, -7.710e22)))

lines_collection = []


for ((prefix, _, label, name), color) in zip(simulations, colors)
    filename = "$(name)_ape_with_time.jld2"
    println(filename)
    if isfile(filename)
        data = load(filename)
        ape = data["APE"]
        time_years = data["time_years"]

        # --- Shift 1/8 interpolated time series
        if name == "eighth_interp"
            time_years .+= last_time_quarter
        end
        

         # --- Shift 1/16 interpolated time series
         if name == "sixteenth_interp"
            time_years .+= (last_time_eighth + last_time_quarter)

        end


        line = lines!(ax, time_years, ape, label=label, linewidth=2.5, color = color)
        push!(lines_collection, line)
    else
        @warn "File not found: $filename"
    end
end

#### code for shifting and zooming in on last 500 years
#=
for ((prefix, _, label, name), color) in zip(simulations, colors)
    filename = "$(name)_ape_with_time.jld2"
    println("Loading $filename")
    if isfile(filename)
        data = load(filename)
        ape = data["APE"]
        time_years = data["time_years"]

        # --- Shifting logic
        if name == "eighth_interp"
            time_years .+= last_time_quarter
        elseif name == "sixteenth_interp"
            time_years .+= last_time_quarter + last_time_eighth
        end

        # --- Trim to last 500 years shown
        total_start_time = last(time_years) - 500
        keep = time_years .>= total_start_time
        time_plot = time_years[keep]
        ape_plot = ape[keep]

        lines!(ax, time_plot, ape_plot, label=label, linewidth=2.5, color=color)
    else
        @warn "File not found: $filename"
    end
end

# --- Add vertical dashed lines at last 1/4° time and last 1/8° interp shifted time
#vlines!(ax, [last_time_quarter], color = :black, linestyle = :dash, linewidth = 2)
#vlines!(ax, [last_time_quarter + last_time_eighth], color = :black, linestyle = :dash, linewidth = 2)
=#
axislegend(ax, labelsize = 20, position = :rc)
display(fig)
save("figures/ape_vs_time_all_shifted.png", fig)
