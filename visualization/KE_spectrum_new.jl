using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using WenoNeverworld.Diagnostics: Spectrum, average_spectra
using CairoMakie
using FFTW
using CUDA

CUDA.device!(1)

# Define function to compute kinetic energy spectrum

function kinetic_energy_spectrum(fields, xrange, yrange)
    # Compute kinetic energy field
    
    Uspec = average_spectra(fields[:u], xrange, yrange)
    Vspec = average_spectra(fields[:v], xrange, yrange)
    
    return Uspec + Vspec
end

# Load data for different grid resolutions
prefixes = ["weno_half_ch", "weno_quarter__ch", "weno_eighth_ch", "weno_sixteen_ch", "weno_thirtytwo_comp"]
dirs = ["/storage4/WenoNeverworldData/half_degree_new/", "/storage4/WenoNeverworldData/quarter_degree_new/", "/storage4/WenoNeverworldData/eighth_degree_new/", "/storage4/WenoNeverworldData/sixteenth_degree_new/", "/storage4/WenoNeverworldData/"]
resolutions = ["1/2", "1/4", "1/8", "1/16", "1/32"]
colors = [:red3, :darkorange, :green, :purple, :navy]


xrange = [11 : 60 * i - 10  for i in [2, 4, 8, 16, 32]]
y50S   = [42  * i : 48  * i for i in [1, 2, 4, 8, 16]]  # channel           #index_50S = closest_lat_index(lats, -50) = 21 ---> 21*2 = 42
y12N   = [160 * i : 170 * i for i in [1, 2, 4, 8, 16]]  #equator            #gives 83 --> 83*2 = 166
y37N   = [210 * i : 220 * i for i in [1, 2, 4, 8, 16]]  #subpolar gyre      #gives 108
y60N   = [260 * i : 270 * i for i in [1, 2, 4, 8, 16]]  #subtropical gyre   


KSPEC2  = Spectrum[] 
KSPEC4  = Spectrum[] 
KSPEC8  = Spectrum[] 
KSPEC16 = Spectrum[] 
#KSPEC32 = Spectrum[] 

i  = 1 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y12N[i], y37N[i], y60N[i]) 
  push!(KSPEC2, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 2 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y12N[i], y37N[i], y60N[i])
  push!(KSPEC4, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

i  = 3 
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 30)
for yr in (y50S[i], y12N[i], y37N[i], y60N[i]) 
  push!(KSPEC8, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)

#=
i  = 4
xr = xrange[i]
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
for yr in (y50S[i], y12N[i], y37N[i], y60N[i]) 
  push!(KSPEC16, kinetic_energy_spectrum(fields, xr, yr))
end

GC.gc(true)
=#
using JLD2

# file_pth = '/home/lcbrock/repository/WenoNeverworld.jl/visualization/weno_sixteen_spectra.jld2'
KSPEC16 = jldopen("weno_sixteen_spectra.jld2", "r") do file
    file["spectra"]
end

fig = Figure(resolution = (2000, 2000), fontsize = 25)

ax1 = Axis(fig[1, 1], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25,
		      title  = L"37^\circ\text{N}",
          width = 500, height = 500,
		      ylabel = L"\text{KE [m^2 s^{-2}]}")

ax2 = Axis(fig[1, 2], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25,
		      title  = L"12^\circ\text{N}",
          width = 500, height = 500)

ax3 = Axis(fig[2, 1], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25,
		      title  = L"50^\circ\text{S}",
          width = 500, height = 500,
          xlabel = L"k \, [ 1/m]",
		      ylabel = L"\text{KE [m^2 s^{-2}]}")

ax4 = Axis(fig[2, 2], xscale = log10, yscale = log10, yticklabelsize = 25, xticklabelsize = 25,
		      title  = L"60^\circ\text{N}",
          width = 500, height = 500,
		      xlabel = L"k \, [ 1/m]")        

grid     = NeverworldGrid(1)
delta50S = grid. Δxᶜᶠᵃ[21]
delta15N = grid. Δxᶜᶠᵃ[82]
delta50N = grid. Δxᶜᶠᵃ[108]
delta60N = grid. Δxᶜᶠᵃ[131]

 F2_50S =  KSPEC2[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F4_50S =  KSPEC4[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F8_50S =  KSPEC8[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F16_50S = KSPEC16[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
#F32_50S = KSPEC32[1].freq[2:end] .* 1 ./ delta50S ./ 3.141592 
 F2_12N =  KSPEC2[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F4_12N =  KSPEC4[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F8_12N =  KSPEC8[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F16_12N = KSPEC16[2].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
#F32_15N = KSPEC32[1].freq[2:end] .* 1 ./ delta15N ./ 3.141592 
 F2_37N =  KSPEC2[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F4_37N =  KSPEC4[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F8_37N =  KSPEC8[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
 F16_37N = KSPEC16[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
#F32_50N = KSPEC32[3].freq[2:end] .* 1 ./ delta50N ./ 3.141592 
F2_60N =  KSPEC2[4].freq[2:end] .* 1 ./ delta60N ./ 3.141592 
F4_60N =  KSPEC4[4].freq[2:end] .* 1 ./ delta60N ./ 3.141592 
F8_60N =  KSPEC8[4].freq[2:end] .* 1 ./ delta60N ./ 3.141592 
F16_60N =  KSPEC16[4].freq[2:end] .* 1 ./ delta60N ./ 3.141592 

lines!(ax1,  F2_37N,  KSPEC2[3].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax1,  F4_37N,  KSPEC4[3].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax1,  F8_37N,  KSPEC8[3].spec[2:end], color = colors[3], linewidth = 2)
lines!(ax1, F16_37N, KSPEC16[3].spec[2:end], color = colors[4], linewidth = 2)
#lines!(ax1, F32_50N, KSPEC32[2].spec[2:end], color = colors[5], linewidth = 1.5)

lines!(ax1, F8_37N[20:end-8], F8_37N[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)#, linestyle = :bold)
#lines!(ax3, F16_50N[10:end-18], F16_50N[10:end-18].^(-2) ./ 10^(13.), color = :black, linestyle = :dashdot)

vlines!(ax1, 1 / 40e3, linestyle = :dash, color = :grey)


lines!(ax2,  F2_12N,  KSPEC2[2].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax2,  F4_12N,  KSPEC4[2].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax2,  F8_12N,  KSPEC8[2].spec[2:end], color = colors[3], linewidth = 2)
lines!(ax2, F16_12N, KSPEC16[2].spec[2:end], color = colors[4], linewidth = 2)
#lines!(ax2, F32_15N, KSPEC32[2].spec[2:end], color = colors[5], linewidth = 1.5)

#lines!(ax2, F8_15N[20:end-8], F8_15N[20:end-8].^(-3) ./ 10^(18), color = :red) #, linestyle = :dashdot)
lines!(ax2, F8_12N[20:end-8], F8_12N[20:end-8].^(-3) ./ 10^(18.5), linewidth = 2.5, color = :black)#, linestyle = :bold)

vlines!(ax2, 1 / 100e3, linestyle = :dash, color = :grey) 

l1 =lines!(ax3,  F2_50S,  KSPEC2[1].spec[2:end], color = colors[1], linewidth = 2, label = L"1/2-\text{degree resolution}")
l2 = lines!(ax3,  F4_50S,  KSPEC4[1].spec[2:end], color = colors[2], linewidth = 2, label = L"1/4-\text{degree resolution}")
l3 = lines!(ax3,  F8_50S,  KSPEC8[1].spec[2:end], color = colors[3], linewidth = 2, label = L"1/8-\text{degree resolution}")
l4 = lines!(ax3, F16_50S, KSPEC16[1].spec[2:end], color = colors[4], linewidth = 2, label = L"1/16-\text{degree resolution}")
#l5 = lines!(ax3, F32_50S, KSPEC32[1].spec[2:end], color = colors[5], linewidth = 1.5, label = L"1/32-\text{degree resolution}")
l6 = lines!(ax3, F8_50S[20:end-8], F8_50S[20:end-8].^(-3) ./ 10^(17.5), linewidth = 2.5, color = :black)

l7 = vlines!(ax3, 1 / 30e3, linestyle = :dash, color = :grey)

######## 60N #########
lines!(ax4,  F2_60N,  KSPEC2[2].spec[2:end], color = colors[1], linewidth = 2)
lines!(ax4,  F4_60N,  KSPEC4[2].spec[2:end], color = colors[2], linewidth = 2)
lines!(ax4,  F8_60N,  KSPEC8[2].spec[2:end], color = colors[3], linewidth = 2)
lines!(ax4, F16_60N, KSPEC16[2].spec[2:end], color = colors[4], linewidth = 2)
#lines!(ax4, F32_60N, KSPEC32[2].spec[2:end], color = colors[5], linewidth = 1.5)
lines!(ax4, F8_60N[20:end-8], F8_12N[20:end-8].^(-3) ./ 10^(18.5), linewidth = 2.5, color = :black)

vlines!(ax4, 1 / 25e3, linestyle = :dash, color = :grey) #need to add correct deformation radius here -- I think this is it?


#leg = Legend(fig[1, 4], ax1)
axislegend(ax1,
  [l1, l2, l3, l4, l6, l7],
  ["1/2°", "1/4°", "1/8°", "1/16°", "-3 slope", L"L_D"], labelsize = 25, position = :rt)
axislegend(ax2,
  [l1, l2, l3, l4, l6, l7],
  ["1/2°", "1/4°", "1/8°", "1/16°", "-3 slope", L"L_D"], labelsize = 25, position = :rt)
axislegend(ax3,
[l1, l2, l3, l4, l6, l7],
["1/2°", "1/4°", "1/8°", "1/16°", "-3 slope", L"L_D"], labelsize = 25, position = :rt)
axislegend(ax4,
  [l1, l2, l3, l4, l6, l7],
  ["1/2°", "1/4°", "1/8°", "1/16°", "-3 slope", L"L_D"], labelsize = 25, position = :rt)

resize_to_layout!(fig)      
display(fig)

CairoMakie.activate!()
CairoMakie.save("figures/ke_spectrum5.png", fig, px_per_unit = 3)
# CairoMakie.save("spectra.eps", fig)

