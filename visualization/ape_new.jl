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

# function to load/save APE and time
using JLD2


function load_ape_with_time(prefix, dir, name; stride=20)
    println("Loading data for: $prefix", "$dir")
    fields = all_fieldtimeseries(prefix, dir; variables=("b",), checkpointer=true)
    time_years = Float64[]
    # Extract and sort the checkpoint keys

    #checkpoints = sort(collect(keys(fields[:b])))
    checkpoints = sort(vec(collect(keys(fields[:b]))))
    ape = Diagnostics.integral_available_potential_energy(fields[:b]; stride)
    
    end_time = length(fields[:b])
    
    for i in 1:stride:end_time
        checkpoint_time = fields[:t][i]
        simulation_time_seconds = checkpoint_time
        simulation_days = simulation_time_seconds / 86400
        sim_years = simulation_days / 365.25

        push!(time_years, sim_years)

        println("Time: $(round(sim_years, digits=3)) years")
    end

      

    # Define the filename for saving
    filename = "$(name)_ape_with_time.jld2"
    println("Saving to: $filename")

    # Save APE and times to a JLD2 file
    jldopen(filename, "w") do file
        write(file, "APE", ape)
        write(file, "time_years", time_years)
    end

    return ape, time_years 
end





simulations = [
   #("weno_half_ch", "/storage4/WenoNeverworldData/half_degree_new/", "1/2∘", "weno_half_ch"),
   # ("weno_quarter__ch", "/storage4/WenoNeverworldData/quarter_degree_new/", "1/4∘", "weno_quarter__ch"),
    #("weno_eighth_checkpoint_iteration2", "/storage4/WenoNeverworldData/eighth_degree_new/", "1/8∘", "eighth2"),
    #("weno_eighth_checkpoint_iteration", "/storage4/WenoNeverworldData/eighth_degree_interp/", "1/8∘ interp", "eighth_interp"),
    ("weno_sixteen_checkpoint_iteration", "/storage4/WenoNeverworldData/sixteenth_degree/", "1/16∘", "sixteenth_interp")
]

for (prefix, dir, _, name) in simulations
    println("Processing simulation: $prefix")
    ape, time_years = load_ape_with_time(prefix, dir, name; stride=20)
    println("Saved: $(name)_ape_with_time.jld2")
end

#=
#plotting starts here

fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
          xlabel = "Time [years]",
          ylabel = L"APE \; (m^5\,s^{-2})",
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
        

        line = lines!(ax, time_years, ape, label=label)
        push!(lines_collection, line)
    else
        @warn "File not found: $filename"
    end
end

axislegend(ax, position = :rb)
display(fig)
save("figures/ape_vs_time_all.png", fig)
=#