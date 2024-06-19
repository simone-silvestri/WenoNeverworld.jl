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
prefixes = ["weno_half_ch", "weno_fourth_ch", "weno_eighth_ch", "weno_sixteenth_ch", "weno_thirtytwo_comp"]
dirs = ["/storage2/WenoNeverworldData/", "/storage2/WenoNeverworldData/", "/storage3/WenoNeverworldData/", "/storage3/WenoNeverworldData/", "/storage4/WenoNeverworldData/"]
resolutions = ["1/2", "1/4", "1/8", "1/16", "1/32"]
colors = [:red3, :darkorange, :green, :purple, :navy]


xrange = [11 : 60 * i - 10  for i in [2, 4, 8, 16, 32]]
y50S   = [42  * i : 48  * i for i in [1, 2, 4, 8, 16]]
y15N   = [160 * i : 170 * i for i in [1, 2, 4, 8, 16]]
y50N   = [210 * i : 220 * i for i in [1, 2, 4, 8, 16]]

KSPEC2  = Spectrum[] 
KSPEC4  = Spectrum[] 
KSPEC8  = Spectrum[] 
KSPEC16 = Spectrum[] 
KSPEC32 = Spectrum[] 

i  = 1 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC2, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 2 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC4, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 3 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC8, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 4
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC16, kinetic_energy_spectrum(fields, xr, yr))
end

i  = 5
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 2)
for yr in (y50S[i], y15N[i], y50N[i]) 
  push!(KSPEC32, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)


fig = Figure(resolution = (2000, 500))

ax1 = Axis(fig[1, 1], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25) 
		      # title  = L"50^\text{o}\text{S}") 
		      # xlabel = L"\text{Wavenumber}", 
		      # ylabel = L"\text{Kinetic Energy spectra}")

ax2 = Axis(fig[1, 2], xscale = log10, yscale = log10,yticklabelsize = 25, xticklabelsize = 25) 
		      # title  = L"12^\text{o}\text{N}") 
		      # xlabel = L"\text{Wavenumber}")

ax3 = Axis(fig[1, 3], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25)
		      # title  = L"37^\text{o}\text{N}")
		      # xlabel = L"\text{Wavenumber}")

grid     = NeverworldGrid(1)
delta50S = grid. Δxᶜᶠᵃ[21]
delta15N = grid. Δxᶜᶠᵃ[82]
delta50N = grid. Δxᶜᶠᵃ[108]

 F2_50S =  KSPEC2[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F4_50S =  KSPEC4[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F8_50S =  KSPEC8[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
F16_50S = KSPEC16[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
F32_50S = KSPEC32[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F2_15N =  KSPEC2[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F4_15N =  KSPEC4[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F8_15N =  KSPEC8[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
F16_15N = KSPEC16[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
F32_15N = KSPEC32[1].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F2_50N =  KSPEC2[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F4_50N =  KSPEC4[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F8_50N =  KSPEC8[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
F16_50N = KSPEC16[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
F32_50N = KSPEC32[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 

lines!(ax1,  F2_50N,  KSPEC2[3].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax1,  F4_50N,  KSPEC4[3].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax1,  F8_50N,  KSPEC8[3].spec[2:end], color = colors[3], linewidth = 2)
lines!(ax1, F16_50N, KSPEC16[3].spec[2:end], color = colors[4], linewidth = 2)
#lines!(ax1, F32_50N, KSPEC32[2].spec[2:end], color = colors[5], linewidth = 1.5)

lines!(ax1, F16_50N[20:end-8], F16_50N[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :bold)
#lines!(ax3, F16_50N[10:end-18], F16_50N[10:end-18].^(-2) ./ 10^(13.), color = :black, linestyle = :dashdot)

vlines!(ax1, 1 / 40e3, linestyle = :dash, color = :grey)


lines!(ax2,  F2_15N,  KSPEC2[2].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax2,  F4_15N,  KSPEC4[2].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax2,  F8_15N,  KSPEC8[2].spec[2:end], color = colors[3], linewidth = 2)
lines!(ax2, F16_15N, KSPEC16[2].spec[2:end], color = colors[4], linewidth = 2)
#lines!(ax2, F32_15N, KSPEC32[2].spec[2:end], color = colors[5], linewidth = 1.5)

#lines!(ax2, F16_15N[20:end-8], F16_15N[20:end-8].^(-3) ./ 10^(17.5), color = :black) #, linestyle = :dashdot)
lines!(ax2, F16_15N[20:end-8], F16_15N[20:end-8].^(-3) ./ 10^(18.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

vlines!(ax2, 1 / 100e3, linestyle = :dash, color = :grey)

l1 =lines!(ax3,  F2_50S,  KSPEC2[1].spec[2:end], color = colors[1], linewidth = 2, label = L"1/2-\text{degree resolution}")
l2 = lines!(ax3,  F4_50S,  KSPEC4[1].spec[2:end], color = colors[2], linewidth = 2, label = L"1/4-\text{degree resolution}")
l3 = lines!(ax3,  F8_50S,  KSPEC8[1].spec[2:end], color = colors[3], linewidth = 2, label = L"1/8-\text{degree resolution}")
l4 = lines!(ax3, F16_50S, KSPEC16[1].spec[2:end], color = colors[4], linewidth = 2, label = L"1/16-\text{degree resolution}")
#l5 = lines!(ax3, F32_50S, KSPEC32[1].spec[2:end], color = colors[5], linewidth = 1.5, label = L"1/32-\text{degree resolution}")

l6 = lines!(ax3, F16_50S[20:end-8], F16_50S[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :dashdot)

l7 = vlines!(ax3, 1 / 30e3, linestyle = :dash, color = :grey)

#leg = Legend(fig[1, 4], ax1)
leg = Legend(fig[1, 4],
    [l1, l2, l3, l4, l6, l7],
    ["1/2°", "1/4°", "1/8°", "1/16°", "-3 slope", ""], labelsize = 25)

display(fig)
#save("ke_spectrum3_latest2.png", fig)

using CairoMakie
CairoMakie.activate!()
CairoMakie.save("ke_spectrum_final2.png", fig, px_per_unit = 3)
# CairoMakie.save("spectra.eps", fig)

