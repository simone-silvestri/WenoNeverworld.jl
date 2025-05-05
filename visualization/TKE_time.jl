using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using CairoMakie
using LaTeXStrings
using JLD2
using CUDA

CUDA.device!(1)

# Define the simulations: (prefix, directory, label)
simulations = [
   ("weno_half_ch", "", "1/2∘", "weno_half_ch"),
   ("weno_quarter__ch", "", "1/4∘", "weno_quarter__ch"),
   ("weno_eighth_checkpoint_iteration", "", "1/8∘", "weno_eighth_ch"),
   ("weno_eighth_checkpoint_iteration", "", "1/8∘ - int", "eighth_interp"),
   ("weno_sixteen_checkpoint_iteration", "", "1/16∘ - int", "sixteenth_interp"),
]

# Define the stride for sampling
stride = 10
#=
# Function to load and compute TKE and corresponding time in years
function load_tke_with_time(prefix, dir; stride=10)
    println("Loading data for: $prefix in $dir")
    fields = all_fieldtimeseries(prefix, dir; variables=("u", "v", "b"), checkpointer=true)
    time_years = Float64[]
    
    # Compute TKE using Oceananigans diagnostics
    tke = Diagnostics.integral_kinetic_energy(fields[:u], fields[:v]; stride=stride)
    checkpoints = sort(vec(collect(keys(fields[:b]))))
    
    end_time = length(fields[:b])
    
    for i in 1:stride:end_time
        checkpoint_time = fields[:t][i]
        simulation_time_seconds = checkpoint_time
        simulation_days = simulation_time_seconds / 86400
        sim_years = simulation_days / 365.25

        push!(time_years, sim_years)

        println("Time: $(round(sim_years, digits=3)) years")
    end
    
    # Save TKE and time to a JLD2 file
    filename = "$(prefix)_tke_with_time.jld2"
    jldopen(filename, "w") do file
        write(file, "TKE", tke)
        write(file, "time_years", time_years)
    end
    
    return tke, time_years
end

# Process each simulation
for (prefix, dir, _) in simulations
    println("Processing simulation: $prefix")
    tke, time_years = load_tke_with_time(prefix, dir; stride=stride)
    println("Saved: $(prefix)_tke_with_time.jld2")
end
=#

# Plotting
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
          xlabel = L"Time [yrs]",
          ylabel = L"KE \; (m^5\,s^{-2})",
          xlabelsize = 25,
          ylabelsize = 25,
          xticklabelsize = 20,
          yticklabelsize = 20)


# --- First, load last time from quarter-degree simulation
quarter_filename = "weno_quarter__ch_tke_with_time.jld2"
quarter_data = load(quarter_filename)
last_time_quarter = last(quarter_data["time_years"])
println("Last 1/4° time (years): ", last_time_quarter)

# --- Now, load last time from eighth-degree simulation
eighth_filename = "eighth_interp_tke_with_time.jld2"
eighth_data = load(eighth_filename)
last_time_eighth = last(eighth_data["time_years"])
println("Last 1/8° time (years): ", last_time_eighth)

colors = [:red3, :darkorange, :green, :navy, :purple]
lines_collection = []

for ((prefix, _, label, name), color) in zip(simulations, colors)
    filename = "$(name)_tke_with_time.jld2"
    println(filename)
    if isfile(filename)
        data = load(filename)
        tke = data["TKE"]
        time_years = data["time_years"]

        # --- Shift 1/8 interpolated time series
        if name == "eighth_interp"
            time_years .+= last_time_quarter
        end
        

         # --- Shift 1/16 interpolated time series
         if name == "sixteenth_interp"
            time_years .+= (last_time_eighth + last_time_quarter)

        end


        line = lines!(ax, time_years, tke, label=label, linewidth=2.5, color = color)
        push!(lines_collection, line)
    else
        @warn "File not found: $filename"
    end
end




# Add legend
axislegend(ax, position = :rb, labelsize = 20)

# Display and save the figure
display(fig)
save("figures/tke_vs_time_all.png", fig)
