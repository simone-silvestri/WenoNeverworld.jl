using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Grids: AbstractGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2
using CUDA

CUDA.device!(1)

variables = ("b")
stride = 1

# Function to load/save APE and time
using JLD2

function load_tke_with_time(prefix, dir, name; stride=5, batch=2)
    println("Loading data for: $prefix, $dir")

    tke_array = Vector{Float64}()
	time_array = Vector{Float64}()
    # Load the time series data using all_fieldtimeseries (this will avoid loading all files at once)
    myfiles, _ = all_filestimeseries(prefix, dir)
    n_files = length(myfiles) # ensure this is divisible by batch
    n_batch = trunc(Int, n_files/batch)
    for n in 1:batch
        start_idx = (n-1)*n_batch + 1
        end_idx = n*n_batch
        fields = all_fieldtimeseries_stride(prefix, dir; variables=("u", "v"), checkpointer=true,start_file=start_idx,end_file=end_idx,stride=stride)

        # Get the total number of time steps
        #end_time = length(myfiles)
        end_time = length(fields[:t])
        println(end_time)

        # Process each checkpoint one by one to avoid memory overload
        for i in 1:end_time
            checkpoint_time = fields[:t][i]
            simulation_time_seconds = checkpoint_time
            simulation_days = simulation_time_seconds / 86400
            sim_years = simulation_days / 365.25

            #ape = Diagnostics.integral_available_potential_energy(fields[:b]; start_time=i,stride=stride,end_time=i)[1]
            tke = Diagnostics.integral_kinetic_energy(fields[:u], fields[:v]; start_time=i,stride=stride,end_time=i)[1]

            println(sim_years, " ", tke)

            # Save the APE and time for this checkpoint immediately to avoid memory issues
            push!(tke_array, tke)
            push!(time_array, sim_years)

            #println(i, " ", size(sim_years), " ", size(ape))
            #println("[$i / $(end_time)] Time = $(round(sim_years, digits=3)) years, APE = $(round(ape, sigdigits=4))")

            # Free memory after each checkpoint to reduce memory usage
            GC.gc()
        end
    end
    println(tke_array)
    println(time_array)
    
    # Open the JLD2 file for saving APE and time data
    filename = "$(name)_tke_with_time.jld2"
    jldopen(filename, "w") do file
	    write(file, "TKE", tke_array)
	    write(file, "time_years", time_array)
    end
end

# Simulation configuration
simulations = [
   #("weno_half_checkpoint_iteration", "/storage4/WenoNeverworldData/half_degree_new/", "1/2∘", "weno_half_ch"), 
   #("weno_quarter__checkpoint_iteration", "/storage4/WenoNeverworldData/quarter_degree_new/", "1/4∘", "weno_quarter__ch"),
   ("weno_eighth_checkpoint_iteration", "/storage4/WenoNeverworldData/eighth_degree_new/", "1/8∘", "weno_eighth_ch"),
   #("weno_eighth_checkpoint_iteration", "/storage4/WenoNeverworldData/eighth_degree_interp/", "1/8∘ int", "eighth_interp"),
   #("weno_sixteen_checkpoint_iteration", "/storage4/WenoNeverworldData/sixteenth_degree/", "1/16∘", "sixteenth_interp")
]

# Loop over each simulation and process it
for (prefix, dir, _, name) in simulations
    println("Processing simulation: $prefix")
    load_tke_with_time(prefix, dir, name; stride=stride,batch=5)
    println("Saved: $(name)_tke_with_time.jld2")
end
