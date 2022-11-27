using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Operators
using FFTW

using Statistics: mean

# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(center)
# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(weno)

# WenoNeverworld.Diagnostics.calculate_fluctuations!(weno,   [:E, :ζ])
# WenoNeverworld.Diagnostics.calculate_fluctuations!(center, [:E, :ζ])

# Calculate energy (Ê) and enstrophy (Ω̂) spectra on a transect in the channel
# using WenoNeverworld.Diagnostics: average_spectra, hann_window

# Êcenter = average_spectra(center[:Efluc], Colon(), 75:75; k = 68)
# Êweno   = average_spectra(weno[:Efluc],   Colon(), 75:75; k = 68)

# Ω̂center = average_spectra(center[:Efluc], Colon(), 215:215; k = 68, windowing = hann_window)
# Ω̂weno   = average_spectra(weno[:Efluc],   Colon(), 215:215; k = 68, windowing = hann_window)

# NxE = length(xnodes(weno[:Efluc][1]))
# Nxζ = length(xnodes(weno[:ζfluc][1]))

# fx = Êcenter.freq[2:end]
# fxE = fftfreq(NxE, 1 / Δxᶜᶜᶜ(1, 1, 1, weno[:Efluc].grid))[1:Int(NxE ÷ 2)] .* 1e3 # in 1/km
# fxζ = fftfreq(Nxζ, 1 / Δxᶠᶠᶜ(1, 1, 1, weno[:Efluc].grid))[1:Int(Nxζ ÷ 2)] .* 1e3 # in 1/km

# fig = Figure()
# ax = Axis(fig[1, 1], title = "energy spectra", yscale = log10, xscale = log10)
# lines!(ax, fxE[2:end], Êcenter.spec[2:end], color = :red,   label = "center momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
# # lines!(ax, fxE[2:end], Êweno.spec[2:end] ./ Êweno.spec[2] .* Êcenter.spec[2], color = :green, linestyle = :dashdot)
# lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
# lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

# ax = Axis(fig[1, 2], title = "enstrophy spectra", yscale = log10, xscale = log10)
# lines!(ax, fxE[2:end], Ω̂center.spec[2:end], color = :red,   label = "center momentum scheme")
# lines!(ax, fxE[2:end], Ω̂weno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
# # lines!(ax, fxE[2:end], Êweno.spec[2:end] ./ Êweno.spec[2] .* Êcenter.spec[2], color = :green, linestyle = :dashdot)
# lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
# lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

# display(fig)

"""
Instantaneous variables
"""

center = all_fieldtimeseries("files_four_centered_new_bathy/neverworld_quarter_centered_snapshots.jld2")
center = limit_timeseries!(center, center[:b].times[end-30:end])

weno = all_fieldtimeseries("files_four_new_bathy/neverworld_quarter_snapshots.jld2")
weno = limit_timeseries!(weno, weno[:b].times[end-30:end])

"""
Time-averaged variables
"""

variables = ("u", "v", "w", "u2", "v2", "ζ", "ζ2")

center_avg = all_fieldtimeseries("files_four_centered_new_bathy/neverworld_quarter_centered_averages.jld2"; variables)
center_avg = limit_timeseries!(center_avg, center_avg[:b].times[end-30:end])

weno_avg = all_fieldtimeseries("files_four_new_bathy/neverworld_quarter_averages.jld2"; variables)
weno_avg = limit_timeseries!(weno_avg, weno_avg[:b].times[end-30:end])

"""
Mean Kinetic energy and Enstrophy
"""

ūweno   = time_average(weno_avg[:u])
ūcenter = time_average(center_avg[:u])
v̄weno   = time_average(weno_avg[:v])
v̄center = time_average(center_avg[:v])
ζ̄weno   = time_average(weno_avg[:ζ])
ζ̄center = time_average(center_avg[:ζ])

WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(center_avg)
WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(weno_avg)

Etotweno   = time_average(weno_avg[:E])
Etotcenter = time_average(center_avg[:E])
Ωtotweno   = time_average(weno_avg[:ζ2])
Ωtotcenter = time_average(center_avg[:ζ2])

Ēweno   = WenoNeverworld.Diagnostics.KineticEnergyField((u = ūweno,   v = v̄weno))
Ēcenter = WenoNeverworld.Diagnostics.KineticEnergyField((u = ūcenter, v = v̄center))
Ω̄weno   = compute!(Field(ζ̄weno^2))
Ω̄center = compute!(Field(ζ̄center^2))

Eflucweno   = FieldTimeSeries{Center, Center, Center}(weno[:u].grid, weno[:u].times)
Efluccenter = FieldTimeSeries{Center, Center, Center}(weno[:u].grid, weno[:u].times)

Ωflucweno   = FieldTimeSeries{Face, Face, Center}(weno[:u].grid, weno[:u].times)
Ωfluccenter = FieldTimeSeries{Face, Face, Center}(weno[:u].grid, weno[:u].times)

for i in 1:length(weno[:u].times)
    wenovel   = (u = compute!(Field(weno[:u][i]   - ūweno)),   v = compute!(Field(weno[:v][i]   - v̄weno)))
    centervel = (u = compute!(Field(center[:u][i] - ūcenter)), v = compute!(Field(center[:v][i] - v̄center)))

    wenovort   = compute!(Field(VerticalVorticityField(weno, i)   - ζ̄weno))  
    centervort = compute!(Field(VerticalVorticityField(center, i) - ζ̄center))

    set!(Eflucweno[i],   WenoNeverworld.Diagnostics.KineticEnergyField(wenovel))
    set!(Efluccenter[i], WenoNeverworld.Diagnostics.KineticEnergyField(centervel))

    set!(Ωflucweno[i],   compute!(Field(wenovort^2)))
    set!(Ωfluccenter[i], compute!(Field(centervort^2)))
end

Eflucmeanweno1   = time_average(Eflucweno)
Eflucmeancenter1 = time_average(Efluccenter)

Ωflucmeanweno1   = time_average(Ωflucweno)
Ωflucmeancenter1 = time_average(Ωfluccenter)

Eflucmeanweno2   = compute!(Field(Etotweno   - Ēweno))
Eflucmeancenter2 = compute!(Field(Etotcenter - Ēcenter))

Ωflucmeanweno2   = compute!(Field(Ωtotweno   - Ω̄weno))
Ωflucmeancenter2 = compute!(Field(Ωtotcenter - Ω̄wcenter))


# IIĒweno   = compute!(Field(Integral(Ēweno,   dims = (1, 3))))
# IIĒcenter = compute!(Field(Integral(Ēcenter, dims = (1, 3))))


# ϕ = ynodes(IIEweno)

# Eflucmeanweno   = WenoNeverworld.Diagnostics.time_average(Eflucweno)
# Eflucmeancenter = WenoNeverworld.Diagnostics.time_average(Efluccenter)

# IIEflucweno   = compute!(Field(Integral(Eflucmeanweno,   dims = (1, 3))))
# IIEfluccenter = compute!(Field(Integral(Eflucmeancenter, dims = (1, 3))))

# IIEtotweno   = compute!(Field(Integral(Etotweno,   dims = (1, 3))))
# IIEtotcenter = compute!(Field(Integral(Etotcenter, dims = (1, 3))))

# fig = Figure()
# ax = Axis(fig[1, 1])

# lines!(ax, ϕ, interior(IIEflucweno  , 1, :, 1), color = :blue)
# lines!(ax, ϕ, interior(IIEfluccenter, 1, :, 1), color = :red)

# ax = Axis(fig[1, 2])

# lines!(ax, ϕ, interior(IIĒweno  , 1, :, 1), color = :blue)
# lines!(ax, ϕ, interior(IIĒcenter, 1, :, 1), color = :red)

# ax = Axis(fig[1, 3])

# lines!(ax, ϕ, interior(IIEtotweno  , 1, :, 1), color = :blue)
# lines!(ax, ϕ, interior(IIEtotcenter, 1, :, 1), color = :red)

# display(fig)

# IEflucweno   = compute!(Field(Integral(Eflucmeanweno,   dims = 3)))
# IEfluccenter = compute!(Field(Integral(Eflucmeancenter, dims = 3)))

# fig = Figure()
# ax  = Axis(fig[1, 1], title = "Depth integreted kinetic energy, WENO")
# heatmap!(ax, log.(interior(IEflucweno,   :, :, 1) .+ 1e-20), colorrange = (0, log(1e3)), colormap = :magma)

# ax  = Axis(fig[1, 2], title = "Depth integreted kinetic energy, Centered")
# heatmap!(ax, log.(interior(IEfluccenter, :, :, 1) .+ 1e-20), colorrange = (0, log(1e3)), colormap = :magma)

# display(fig)

# using WenoNeverworld.Diagnostics: calculate_z★_diagnostics, calculate_Γ²_diagnostics

# zw = calculate_z★_diagnostics(  weno_avg[:b])
# zc = calculate_z★_diagnostics(center_avg[:b])

# Γw = calculate_Γ²_diagnostics(zw,   weno_avg[:b])
# Γc = calculate_Γ²_diagnostics(zc, center_avg[:b])

# εw = propagate_on_fieldtimeseries(zw,   weno_avg[:b], Γw; func = (x, y, z) -> x * y + z, nargs = 3)
# εc = propagate_on_fieldtimeseries(zc, center_avg[:b], Γc; func = (x, y, z) -> x * y + z, nargs = 3)

# using Oceananigans.Grids: znodes, ynodes

# Bw = calc_zonal_mean(weno_avg[:b])
# Bc = calc_zonal_mean(center_avg[:b])

# ϕ = ynodes(weno_avg[:b][1])
# z = znodes(weno_avg[:b][1])

# fig = Figure()
# ax  = Axis(fig[1, 1], title = "Zonally averaged buoyancy")
# contourf!(ax, ϕ, z, Bw[1, :, :], levels = (0:0.002:0.06))
# contour!(ax,  ϕ, z, Bc[1, :, :], levels = (0:0.002:0.06), color = :black)

# display(fig)

