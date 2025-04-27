using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using WenoNeverworld.Diagnostics: Spectrum, average_spectra
using CairoMakie
using JLD2, Statistics, Printf

# Constants
ρ₀ = 1025  # kg/m³
c_p = 3994  # J/(kg·K)

# List of simulation directories and labels
simulations = [
    (dir="/storage4/WenoNeverworldData/half_degree_new/", label="1/2°"),
    (dir="/storage4/WenoNeverworldData/quarter_degree_new/", label="1/4°"),
    (dir="/storage4/WenoNeverworldData/eighth_degree_new/", label="1/8°")
]

# Initialize figure
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1, 1], xlabel="Latitude [°]", ylabel="MHT [PW]")

for (sim_dir, label) in simulations
    # Load simulation data using checkpoint files
    # Get the time series of fields for u, v, and T
    fields = all_fieldtimeseries("weno_half_ch", sim_dir; variables=("v", "b"), checkpointer=true, number_files=25)
    
    # Get grid and coordinates
    grid = fields.grid
    φ = grid.yC  # Latitude is yC

    # Time-averaging over the loaded field timeseries
    v̄ = time_average(fields[:v])
    b̄ = time_average(fields[:b])

    # Compute MHT as a function of latitude
    MHT = zeros(length(φ))
    for j in 1:length(φ)
        # Extract v and T at latitude j
        v_slice = interior(v̄, :, j, :)
        b_slice = interior(b̄, :, j, :)

        # Compute integrand v * T
        integrand = v_slice .* b_slice

        # Integrate over the x and z dimensions
        dx = xxxxxx
        dz = Δz
        integral = sum(integrand .* dx .* dz)

        # Multiply by constants and convert to petawatts (PW)
        MHT[j] = (ρ₀ * c_p * integral) / 1e15  # Convert to PW
    end

    # Plot MHT vs latitude
    lines!(ax, φ, MHT, label=label)
end

axislegend(ax)
save("figures/MHT_comparison.png", fig)
