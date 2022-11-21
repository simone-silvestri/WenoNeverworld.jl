using WenoNeverworld
using FFTW

using WenoNeverworld.Diagnostics

# centered = all_fieldtimeseries("files_four_centered/neverworld_quarter_centered_snapshots.jld2")
# centered = limit_timeseries!(centered, centered[:b].times[end-20:end])

# weno = all_fieldtimeseries("files_four/neverworld_quarter_snapshots.jld2")
# weno = limit_timeseries!(weno, weno[:b].times[end-20:end])

# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(centered)
# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(weno)

# using WenoNeverworld.Diagnostics: average_spectra

# Êcenter = average_spectra(centered[:E], Colon(), 78:82; k = 65)
# Êweno   = average_spectra(weno[:E], Colon(), 78:82; k = 65)

# Ω̂center = average_spectra(centered[:ζ], Colon(), 78:82; k = 65)
# Ω̂weno   = average_spectra(weno[:ζ], Colon(), 78:82; k = 65)

# NxE = length(xnodes(weno[:E][1]))
# Nxζ = length(xnodes(weno[:ζ][1]))

# fx = Êcenter.freq[2:end]
# fxE = fftfreq(NxE, 1 / Δxᶜᶜᶜ(1, 1, 1, weno[:E].grid))[1:Int(NxE ÷ 2)] .* 1e3 # in 1/km
# fxζ = fftfreq(Nxζ, 1 / Δxᶠᶠᶜ(1, 1, 1, weno[:E].grid))[1:Int(Nxζ ÷ 2)] .* 1e3 # in 1/km

# fig = Figure()
# ax = Axis(fig[1, 1], title = "energy spectra", yscale = log10, xscale = log10)
# lines!(ax, fxE[2:end], Êcenter.spec[2:end], color = :red, label = "Centered momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end], color = :blue, label = "WENO momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end] ./ Êweno.spec[2] .* Êcenter.spec[2], color = :green, linestyle = :dashdot)
# lines!(ax, fxE[10:end-51],   0.001 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
# lines!(ax, fxE[end-100:end], 0.007 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

# ax = Axis(fig[1, 2], title = "enstrophy spectra", yscale = log10, xscale = log10)
# lines!(ax, fxζ[2:end], Ω̂center.spec[2:end], color = :red,  label = "Centered momentum scheme")
# lines!(ax, fxζ[2:end], Ω̂weno.spec[2:end],   color = :blue, label = "WENO momentum scheme")
# lines!(ax, fxζ[2:end], Ω̂weno.spec[2:end] ./ Ω̂weno.spec[2] .* Ω̂center.spec[2], color = :green, linestyle = :dashdot)
# lines!(ax, fxζ[end-110:end-20], 8e-12 .* fx[end-110:end-20].^(-1), color = :grey, linewidth = 2, linestyle = :dashdot)

# display(fig)

center_avg = all_fieldtimeseries("files_four_centered/neverworld_quarter_centered_averages.jld2")
center_avg = limit_timeseries!(centered, centered[:b].times[end-20:end])

weno_avg = all_fieldtimeseries("files_four/neverworld_quarter_averages.jld2")
weno_avg = limit_timeseries!(weno, weno[:b].times[end-20:end])

using WenoNeverworld.Diagnostics: calculate_z★_diagnostics, calculate_Γ²_diagnostics

zw = calculate_z★_diagnostics(  weno_avg[:b])
zc = calculate_z★_diagnostics(center_avg[:b])

Γw = calculate_Γ²_diagnostics(zw,   weno_avg[:b])
Γc = calculate_Γ²_diagnostics(zc, center_avg[:b])

εw = propagate_on_fieldtimeseries(zw,   weno_avg[:b], Γw; func = (x, y, z) -> x * y + z, nargs = 3)
εc = propagate_on_fieldtimeseries(zc, center_avg[:b], Γc; func = (x, y, z) -> x * y + z, nargs = 3)