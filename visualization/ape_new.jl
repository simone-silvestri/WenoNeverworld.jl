using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2
using CUDA

CUDA.device!(2)

function simulated_time_from_checkpoint(filepath)
    jldopen(filepath, "r") do file
        clock = file["clock"]
        time_seconds = clock.time
        return time_seconds / 86400 / 365.25  # Convert to years
    end
end

function load_ape_with_time(prefix, dir; stride=20)
    println("Loading data for: $prefix")
    fields = all_fieldtimeseries(prefix, dir; variables=("b",), checkpointer=true)
    checkpoints = sort(vec(collect(keys(fields[:b]))))

    times_years = Float64[]
    ape_values = Float64[]

    for checkpoint in checkpoints[1:stride:end]
        checkpoint_data = fields[:b][checkpoint]
        checkpoint_path = checkpoint_data.filepath

        if checkpoint_path === nothing || !isfile(checkpoint_path)
            @warn "Missing or invalid checkpoint file: $checkpoint_path"
            continue
        end

        sim_years = simulated_time_from_checkpoint(checkpoint_path)
        ape = Diagnostics.integral_available_potential_energy(checkpoint_data)

        push!(times_years, sim_years)
        push!(ape_values, ape)

        println("Checkpoint: $(basename(checkpoint_path)) → Time: $(round(sim_years, digits=3)) years, APE: $ape")
    end

    filename = "$(prefix)_ape_with_time.jld2"
    jldopen(filename, "w") do file
        write(file, "APE", ape_values)
        write(file, "time_years", times_years)
    end

    return ape_values, times_years
end

simulations = [
    ("weno_half_ch", "/storage4/WenoNeverworldData/half_degree_new/", "1/2∘"),
    ("weno_eighth_checkpoint_iteration3", "/storage4/WenoNeverworldData/eighth_degree_new/", "1/8∘"),
    ("weno_quarter__ch", "/storage4/WenoNeverworldData/quarter_degree_new/", "1/4∘")
]

for (prefix, dir, _) in simulations
    println("Processing simulation: $prefix")
    _, _ = load_ape_with_time(prefix, dir; stride=20)
end

fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1],
    xlabel = "Time [years]",
    ylabel = L"APE \; (m^5\,s^{-2})",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 20,
    yticklabelsize = 20
)

for (prefix, _, label) in simulations
    filename = "$(prefix)_ape_with_time.jld2"
    if isfile(filename)
        data = load(filename)
        lines!(ax, data["time_years"], data["APE"], label=label)
    else
        @warn "Data file not found: $filename"
    end
end

axislegend(ax, position=:rb)
save("figures/ape_vs_time_all_simulations.png", fig)
display(fig)
