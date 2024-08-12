using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using WenoNeverworld.Diagnostics: Spectrum, average_spectra
using GLMakie
using FFTW

# Define function to compute kinetic energy spectrum

function kinetic_energy_spectrum(fields, xrange, yrange)
    # Compute kinetic energy field
    
    Uspec = average_spectra(fields[:u], xrange, yrange)
    Vspec = average_spectra(fields[:v], xrange, yrange)
    
    return Uspec + Vspec
end

# Load data for different grid resolutions
prefixes = ["weno_quarter_ch", "weno_quarter_higher_diffusivity_ch", "weno_half_original_ch", "weno_half_larger_diffusivity_ch"]
dirs = ["/storage4/WenoNeverworldData/quarter_degree/", "/storage4/WenoNeverworldData/quarter_degree/", "/storage4/WenoNeverworldData/half_degree/", "/storage4/WenoNeverworldData/half_degree/"]
resolutions = ["1/4 original", "1/4 larger diff", "1/2 original", "1/2 larger diff"]
colors = [:blue, :red, :green, :purple]


xrange = [11 : 60 * i - 10  for i in [4, 4, 2, 2]]  #remember to change these for different resolutions
y50S   = [42  * i : 48  * i for i in [2, 2, 1, 1]]
y15N   = [160 * i : 170 * i for i in [2, 2, 1, 1]]
y50N   = [210 * i : 220 * i for i in [2, 2, 1, 1]]

KSPEC1  = Spectrum[] 
KSPEC2  = Spectrum[] 
KSPEC3  = Spectrum[] 
KSPEC4  = Spectrum[] 

num_files = 10 # could change this so it doesn't return num_files + 1 files
print()

i  = 1 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = num_files)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC1, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 2 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = num_files)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC2, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)


i  = 3 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = num_files)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC3, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 4
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = num_files)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC4, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

fig = Figure(resolution = (2000, 500))

ax1 = Axis(fig[1, 1], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25, 
		    title  = "50 S", xlabel = "k", ylabel = "Kinetic Energy spectra")

ax2 = Axis(fig[1, 2], xscale = log10, yscale = log10,yticklabelsize = 25, xticklabelsize = 25, 
		    title  = "12 N", xlabel = "k")

ax3 = Axis(fig[1, 3], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25,
		      title  = "37 N", xlabel = "k")

grid     = NeverworldGrid(1)
delta50S = grid. Δxᶜᶠᵃ[21]
delta15N = grid. Δxᶜᶠᵃ[82]
delta50N = grid. Δxᶜᶠᵃ[108]

# Function to filter out non-positive values
function filter_non_positive!(x, y)
  mask = x .> 0 .&& y .> 0
  return x[mask], y[mask]
end

F2_50S, KSPEC1_50S = filter_non_positive!(KSPEC1[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592, KSPEC1[1].spec[2:end])
F4_50S, KSPEC2_50S = filter_non_positive!(KSPEC2[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592, KSPEC2[1].spec[2:end])
F8_50S, KSPEC3_50S = filter_non_positive!(KSPEC3[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592, KSPEC3[1].spec[2:end])
F16_50S, KSPEC4_50S = filter_non_positive!(KSPEC4[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592, KSPEC4[1].spec[2:end])

F2_15N, KSPEC1_15N = filter_non_positive!(KSPEC1[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592, KSPEC1[2].spec[2:end])
F4_15N, KSPEC2_15N = filter_non_positive!(KSPEC2[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592, KSPEC2[2].spec[2:end])
F8_15N, KSPEC3_15N = filter_non_positive!(KSPEC3[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592, KSPEC3[2].spec[2:end])
F16_15N, KSPEC4_15N = filter_non_positive!(KSPEC4[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592, KSPEC4[2].spec[2:end])

F2_50N, KSPEC1_50N = filter_non_positive!(KSPEC1[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592, KSPEC1[3].spec[2:end])
F4_50N, KSPEC2_50N = filter_non_positive!(KSPEC2[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592, KSPEC2[3].spec[2:end])
F8_50N, KSPEC3_50N = filter_non_positive!(KSPEC3[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592, KSPEC3[3].spec[2:end])
F16_50N, KSPEC4_50N = filter_non_positive!(KSPEC4[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592, KSPEC4[3].spec[2:end])

lines!(ax1, F2_50N, KSPEC1_50N, color = colors[1], linewidth = 2)
lines!(ax1, F4_50N, KSPEC2_50N, color = colors[2], linewidth = 2)
lines!(ax1, F8_50N, KSPEC3_50N, color = colors[3], linewidth = 2)
lines!(ax1, F16_50N, KSPEC4_50N, color = colors[4], linewidth = 2)
#lines!(ax1, F16_50N[20:end-8], F16_50N[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

#vlines!(ax1, 1 / 100e3, linestyle = :dash, color = :grey)

lines!(ax2, F2_15N, KSPEC1_15N, color = colors[1], linewidth = 2)
lines!(ax2, F4_15N, KSPEC2_15N, color = colors[2], linewidth = 2)
lines!(ax2, F8_15N, KSPEC3_15N, color = colors[3], linewidth = 2)
lines!(ax2, F16_15N, KSPEC4_50N, color = colors[4], linewidth = 2)
#lines!(ax2, F16_15N[20:end-8], F16_15N[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

#vlines!(ax2, 1 / 100e3, linestyle = :dash, color = :grey)

#=

 F2_50S =  KSPEC1[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F4_50S =  KSPEC2[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F8_50S =  KSPEC3[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F16_50S = KSPEC4[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592

 F2_15N =  KSPEC1[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F4_15N =  KSPEC2[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F8_15N =  KSPEC3[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F16_15N = KSPEC4[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 

 F2_50N =  KSPEC1[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F4_50N =  KSPEC2[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F8_50N =  KSPEC3[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F16_50N = KSPEC4[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592

lines!(ax1,  F2_50N,  KSPEC1[3].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax1,  F4_50N,  KSPEC2[3].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax1,  F8_50N,  KSPEC3[3].spec[2:end], color = colors[3], linewidth = 2)
#lines!(ax1, F16_50N,  KSPEC4[3].spec[2:end], color = colors[4], linewidth = 2)

lines!(ax1, F16_50N[20:end-8], F16_50N[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :bold)

vlines!(ax1, 1 / 40e3, linestyle = :dash, color = :grey)

lines!(ax2,  F2_15N,  KSPEC1[2].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax2,  F4_15N,  KSPEC2[2].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax2,  F8_15N,  KSPEC3[2].spec[2:end], color = colors[3], linewidth = 2)
#lines!(ax2, F16_15N,  KSPEC4[2].spec[2:end], color = colors[4], linewidth = 2)

lines!(ax2, F16_15N[20:end-8], F16_15N[20:end-8].^(-3) ./ 10^(18.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

vlines!(ax2, 1 / 100e3, linestyle = :dash, color = :grey)
=#
l1 =lines!(ax3,   F2_50S,  KSPEC1[1].spec[2:end], color = colors[1], linewidth = 2, label = L"1/2-\text{degree resolution}")
l2 = lines!(ax3,  F4_50S,  KSPEC2[1].spec[2:end], color = colors[2], linewidth = 2, label = L"1/4-\text{degree resolution}")
l3 = lines!(ax3,  F8_50S,  KSPEC3[1].spec[2:end], color = colors[3], linewidth = 2, label = L"1/8-\text{degree resolution}")
l4 = lines!(ax3, F16_50S,  KSPEC4[1].spec[2:end], color = colors[4], linewidth = 2, label = L"1/16-\text{degree resolution}")

l5 = lines!(ax3, F4_50S[20:end-8], F4_50S[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

#l6 = vlines!(ax3, 1 / 100e3, linestyle = :dash, color = :grey)

leg = Legend(fig[1, 4],
    [l1, l2, l3, l4, l5],
    ["1/4 original", "1/4 larger diff", "1/2 original", "1/2 larger diff", "-3 slope"], labelsize=25)

display(fig)
#save("ke_spectrum3_latest2.png", fig)

using CairoMakie
CairoMakie.activate!()
CairoMakie.save("ke_spec2.png", fig, px_per_unit = 3)

#ke_spec1 = last 20 files
#ke_spec2 = last 10 files

