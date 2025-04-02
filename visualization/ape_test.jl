using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2
using CUDA

CUDA.device!(1)

variables = ("b")
stride = 20

# Function to get all checkpoint files from a directory with stride
function get_checkpoint_files(dir, stride)
    # List all files in the directory
    all_files = readdir(dir)
    
    # Filter files that contain "checkpoint" in their name (you can modify the pattern if needed)
    checkpoint_files = filter(f -> contains(f, "checkpoint"), all_files)
    
    # Sort the files based on their iteration (extract the number from the filename)
    checkpoint_files = sort(checkpoint_files, by = f -> parse(Int, match(r"\d+", f).match))
    
    # Select files based on the stride (every 20th file)
    selected_files = checkpoint_files[1:stride:end]
    
    # Return full paths of the selected checkpoint files
    return [joinpath(dir, f) for f in selected_files]
end

# Function to extract simulated time from checkpoint file
function simulated_time_from_checkpoint(checkpoint_file)
    file = jldopen(checkpoint_file, "r")
    try
        clock = file["clock"]
        simulation_time_seconds = clock.time
        simulation_days = simulation_time_seconds / 86400
        simulation_years = simulation_days / 365.25
        return simulation_days, simulation_years
    finally
        close(file)
    end
end

# Function to extract APE from a simulation
function extract_ape_from_simulation(prefix_simulation, dir, checkpoint_files, stride)
    ape_data = []
    times = []
    
    for checkpoint_file in checkpoint_files
        fields = all_fieldtimeseries(prefix_simulation, dir; variables, checkpointer = true)
        ape = Diagnostics.integral_available_potential_energy(fields[:b]; stride)
        
        # Extract simulation time from the checkpoint file
        sim_days, sim_years = simulated_time_from_checkpoint(checkpoint_file)
        
        # Store APE data and corresponding simulated time
        push!(ape_data, ape)
        push!(times, sim_days)  # or sim_years if you prefer
        
    end
    
    return ape_data, times
end

# Define directories for each simulation resolution
dir_half = "/storage4/WenoNeverworldData/half_degree_new/"
dir_fourth = "/storage4/WenoNeverworldData/quarter_degree_new/"
dir_eighth = "/storage4/WenoNeverworldData/eighth_degree_new/"

# Automatically get the checkpoint files for each resolution with stride
checkpoint_files_half = get_checkpoint_files(dir_half, stride)
checkpoint_files_fourth = get_checkpoint_files(dir_fourth, stride)
checkpoint_files_eighth = get_checkpoint_files(dir_eighth, stride)

# Extract data for half-degree simulations
prefix_half = "weno_half_ch"
ape_half, time_half = extract_ape_from_simulation(prefix_half, dir_half, checkpoint_files_half, stride)

# Extract data for quarter-degree simulations
prefix_fourth = "weno_quarter_ch"
ape_fourth, time_fourth = extract_ape_from_simulation(prefix_fourth, dir_fourth, checkpoint_files_fourth, stride)

# Extract data for eighth-degree simulations
prefix_eighth = "weno_eighth_ch"
ape_eighth, time_eighth = extract_ape_from_simulation(prefix_eighth, dir_eighth, checkpoint_files_eighth, stride)

# Now plot the data using CairoMakie
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1], xlabel="Simulated Time (days)", ylabel=L"[m^5/s^2]", title="APE vs. Simulated Time")

# Plot APE data for each resolution
lines!(ax, time_half, ape_half, color=:blue, label="1/2 degree")
lines!(ax, time_fourth, ape_fourth, color=:red, label="1/4 degree")
lines!(ax, time_eighth, ape_eighth, color=:green, label="1/8 degree")

# Add a legend
Legend(fig[1, 2], [l1, l2, l3], ["1/2 degree", "1/4 degree", "1/8 degree"], position=:right)

# Display the plot
display(fig)

# Save the figure
save("ape_vs_time_comparison.png", fig)
