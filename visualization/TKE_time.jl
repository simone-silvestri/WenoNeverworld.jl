using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using CairoMakie
using LaTeXStrings
using JLD2
using CUDA

# Set the CUDA device
CUDA.device!(2)

# Define the simulations: (prefix, directory, label)
simulations = [
    ("weno_half_ch", "/storage4/WenoNeverworldData/half_degree_new/", "1/2∘"),
    ("weno_quarter__ch", "/storage4/WenoNeverworldData/quarter_degree_new/", "1/4∘"),
    ("weno_eighth_checkpoint_iteration3", "/storage4/WenoNeverworldData/eighth_degree_new/", "1/8∘")
]

# Define the stride for sampling
stride = 5

# Function to load and compute TKE and corresponding time in years
function load_tke_with_time(prefix, dir; stride=5)
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

# Plotting
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
          xlabel = "Time [years]",
          ylabel = L"TKE \; (m^5\,s^{-2})",
          xlabelsize = 20,
          ylabelsize = 20,
          xticklabelsize = 20,
          yticklabelsize = 20)

# Plot TKE for each simulation
for (prefix, _, label) in simulations
    filename = "$(prefix)_tke_with_time.jld2"
    
    if isfile(filename)
        data = load(filename)
        tke = data["TKE"]
        time_years = data["time_years"]
        
        lines!(ax, time_years, tke, label=label)
    else
        @warn "File not found: $filename"
    end
end

# Add legend
axislegend(ax, position = :rb)

# Display and save the figure
display(fig)
save("figures/tke_vs_time1.png", fig)
