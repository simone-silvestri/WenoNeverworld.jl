using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2

# --- List of simulations
simulations = [
   ("weno_half_ch", "", "1/2∘", "weno_half_ch"),
   ("weno_quarter__ch", "", "1/4∘", "weno_quarter__ch"),
   ("weno_eighth_checkpoint_iteration2", "", "1/8∘", "eighth2"),
   ("weno_eighth_checkpoint_iteration", "", "1/8∘ interp", "eighth_interp"),
   ("weno_sixteen_checkpoint_iteration", "", "1/16∘ interp", "sixteenth_interp"),
]

# --- First, load last time from quarter-degree simulation
quarter_filename = "weno_quarter__ch_ape_with_time.jld2"
quarter_data = load(quarter_filename)
last_time_quarter = last(quarter_data["time_years"])

println("Last 1/4° time (years): ", last_time_quarter)

# --- Now, load last time from eighth-degree simulation
eighth_filename = "eighth_interp"
eighth_data = load(eighth_filename)
last_time_eighth = last(eighth_data["time_years"])

println("Last 1/8° time (years): ", last_time_eighth)

# --- Now plot everything
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
          xlabel = L"Time [years]",
          ylabel = L"APE \; [m^5\,s^{-2}]",
          xlabelsize = 20,
          ylabelsize = 20,
          xticklabelsize = 20,
          yticklabelsize = 20)

lines_collection = []

for (prefix, _, label, name) in simulations
    filename = "$(name)_ape_with_time.jld2"

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
            time_years .+= last_time_eighth
        end


        line = lines!(ax, time_years, ape, label=label)
        push!(lines_collection, line)
    else
        @warn "File not found: $filename"
    end
end



axislegend(ax, position = :rb)
display(fig)
save("figures/ape_vs_time_all_shifted2.png", fig)
