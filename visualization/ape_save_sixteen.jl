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

# Function to load/save APE and time
using JLD2

function load_ape_with_time(prefix, dir, name; stride=5)
    println("Loading data for: $prefix, $dir")

    # Load the time series data using all_fieldtimeseries (this will avoid loading all files at once)
    fields = all_fieldtimeseries(prefix, dir; variables=("b",), checkpointer=true)

    # Open the JLD2 file for saving APE and time data
    filename = "$(name)_ape_with_time.jld2"
    jldopen(filename, "w") do file
        # Get the total number of time steps
        end_time = length(fields[:b])

        # Process each checkpoint one by one to avoid memory overload
        for i in 1:stride:end_time
            try
                # Extract the buoyancy and time for the current checkpoint
                b = fields[:b][i]
                time = fields[:t][i]

                # Calculate APE
                ape = Diagnostics.integral_available_potential_energy(b)

                sim_years = time / 86400 / 365.25

                # Save the APE and time for this checkpoint immediately to avoid memory issues
                write(file, "APE_$i", ape)
                write(file, "time_years_$i", sim_years)

                println("[$i / $(end_time)] Time = $(round(sim_years, digits=3)) years, APE = $(round(ape, sigdigits=4))")

                # Free memory after each checkpoint to reduce memory usage
                GC.gc()

            catch e
                @warn "Skipping checkpoint $i due to error: $e"
                continue
            end
        end
    end
end

# Simulation configuration
simulations = [
   ("weno_sixteen_checkpoint_iteration", "/storage4/WenoNeverworldData/sixteenth_degree/", "1/16∘", "sixteenth_interp")
]

# Loop over each simulation and process it
for (prefix, dir, _, name) in simulations
    println("Processing simulation: $prefix")
    load_ape_with_time(prefix, dir, name; stride=5)
    println("Saved: $(name)_ape_with_time.jld2")
end