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
stride = 5

# Function to load/save APE and time
using JLD2

function load_ape_with_time(prefix, dir, name; stride=5)
    println("Loading data for: $prefix, $dir")

    # Load the time series data using all_fieldtimeseries (this will avoid loading all files at once)
    fields = all_fieldtimeseries(prefix, dir; variables=("b",), checkpointer=true)

    #myfiles, numbers = all_filestimeseries(prefix, dir)
    
    # Get the total number of time steps
    #end_time = length(myfiles)
    end_time = length(fields[:t])
    println(end_time)

    ape_array = Vector{Float64}()
	time_array = Vector{Float64}()
    
    # Process each checkpoint one by one to avoid memory overload
    for i in 1:stride:end_time
        #try
            #=
            # Extract the buoyancy and time for the current checkpoint
            time = jldopen(dir*myfiles[i]*"2")["clock"].time
            println(time)
            grid = jldopen(dir*myfiles[i]*"2")["grid"]
            Hx = grid.underlying_grid.Hx
            Hy = grid.underlying_grid.Hy
            Hz = grid.underlying_grid.Hz
            b = jldopen(dir*myfiles[i]*"2")["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
            vol = VolumeField(grid)

            # Calculate APE
            αe = Diagnostics.compute_ape_density(Field(b))
            ape = sum(compute!(Field(αe * vol)))

            println(size(ape))

            sim_years = time / 86400 / 365.25
            =#
            checkpoint_time = fields[:t][i]
            simulation_time_seconds = checkpoint_time
            simulation_days = simulation_time_seconds / 86400
            sim_years = simulation_days / 365.25

            ape = Diagnostics.integral_available_potential_energy(fields[:b]; start_time=i,stride=stride,end_time=i)[1]

            println(sim_years, " ", ape)

            # Save the APE and time for this checkpoint immediately to avoid memory issues
            push!(ape_array, ape)
            push!(time_array, sim_years)

            #println(i, " ", size(sim_years), " ", size(ape))
            #println("[$i / $(end_time)] Time = $(round(sim_years, digits=3)) years, APE = $(round(ape, sigdigits=4))")

            # Free memory after each checkpoint to reduce memory usage
            #GC.gc()

        #catch e
        #    @warn "Skipping checkpoint $i due to error: $e"
        #    continue
        #end
    end
    println(ape_array)
    println(time_array)
    
    # Open the JLD2 file for saving APE and time data
    filename = "$(name)_ape_with_time.jld2"
    jldopen(filename, "w") do file
	    write(file, "APE", ape_array)
	    write(file, "time_years", time_array)
    end
end

# Simulation configuration
simulations = [
   #("weno_sixteen_checkpoint_iteration", "/storage4/WenoNeverworldData/sixteenth_degree/", "1/16∘", "sixteenth_interp"), 
   ("weno_eighth_checkpoint_iteration", "/storage4/WenoNeverworldData/eighth_degree_new/", "1/8∘", "weno_eighth")
]

# Loop over each simulation and process it
for (prefix, dir, _, name) in simulations
    println("Processing simulation: $prefix")
    load_ape_with_time(prefix, dir, name; stride=stride)
    println("Saved: $(name)_ape_with_time.jld2")
end
