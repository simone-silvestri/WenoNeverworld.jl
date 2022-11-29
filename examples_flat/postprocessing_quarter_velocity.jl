using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Operators
using GLMakie
using FFTW

using Statistics: mean

"""
Instantaneous variables
"""

variables = ("u", "v")

center = all_fieldtimeseries("files_four_centered_new_bathy/neverworld_quarter_centered_snapshots.jld2")
center = limit_timeseries!(center, center[:b].times[end-30:end])

weno = all_fieldtimeseries("files_four_new_bathy/neverworld_quarter_snapshots.jld2")
weno = limit_timeseries!(weno, weno[:b].times[end-30:end])

trac = all_fieldtimeseries("files_four_isopycnal_new_bathy/neverworld_quarter_isopycnal_snapshots.jld2")
trac = limit_timeseries!(trac, trac[:b].times[end-30:end])

"""
Time-averaged variables
"""

variables = ("u", "v", "w", "u2", "v2", "ζ", "ζ2", "b2")

center_avg = all_fieldtimeseries("files_four_centered_new_bathy/neverworld_quarter_centered_averages.jld2"; variables)
center_avg = limit_timeseries!(center_avg, center_avg[:u].times[end-30:end])

weno_avg = all_fieldtimeseries("files_four_new_bathy/neverworld_quarter_averages.jld2"; variables)
weno_avg = limit_timeseries!(weno_avg, weno_avg[:u].times[end-30:end])

trac_avg = all_fieldtimeseries("files_four_isopycnal_new_bathy/neverworld_quarter_isopycnal_averages.jld2"; variables)
trac_avg = limit_timeseries!(trac_avg, trac_avg[:u].times[end-30:end])

"""
Mean Kinetic energy and Enstrophy
"""

ūweno   = time_average(weno_avg[:u])
ūtrac   = time_average(trac_avg[:u])
ūcenter = time_average(center_avg[:u])
v̄weno   = time_average(weno_avg[:v])
v̄trac   = time_average(trac_avg[:v])
v̄center = time_average(center_avg[:v])
ζ̄weno   = time_average(weno_avg[:ζ])
ζ̄trac   = time_average(trac_avg[:ζ])
ζ̄center = time_average(center_avg[:ζ])

WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(center_avg)
WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(weno_avg)
WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(trac_avg)

Etotweno   = time_average(weno_avg[:E])
Etottrac   = time_average(trac_avg[:E])
Etotcenter = time_average(center_avg[:E])
Ωtotweno   = time_average(weno_avg[:ζ2])
Ωtottrac   = time_average(trac_avg[:ζ2])
Ωtotcenter = time_average(center_avg[:ζ2])

Ēweno   = WenoNeverworld.Diagnostics.KineticEnergyField((u = ūweno,   v = v̄weno))
Ētrac   = WenoNeverworld.Diagnostics.KineticEnergyField((u = ūtrac,   v = v̄trac))
Ēcenter = WenoNeverworld.Diagnostics.KineticEnergyField((u = ūcenter, v = v̄center))
Ω̄weno   = compute!(Field(ζ̄weno^2))
Ω̄trac   = compute!(Field(ζ̄trac^2))
Ω̄center = compute!(Field(ζ̄center^2))

Eflucweno   = FieldTimeSeries{Center, Center, Center}(weno[:u].grid, weno[:u].times)
Efluctrac   = FieldTimeSeries{Center, Center, Center}(weno[:u].grid, weno[:u].times)
Efluccenter = FieldTimeSeries{Center, Center, Center}(weno[:u].grid, weno[:u].times)

Ωflucweno   = FieldTimeSeries{Face, Face, Center}(weno[:u].grid, weno[:u].times)
Ωfluctrac   = FieldTimeSeries{Face, Face, Center}(weno[:u].grid, weno[:u].times)
Ωfluccenter = FieldTimeSeries{Face, Face, Center}(weno[:u].grid, weno[:u].times)

using WenoNeverworld.Diagnostics: VerticalVorticityField

for i in 1:length(weno[:u].times)
    @info "doing time $i"
    wenovel   = (u = compute!(Field(weno[:u][i]   - ūweno)),   v = compute!(Field(weno[:v][i]   - v̄weno)))
    tracvel   = (u = compute!(Field(trac[:u][i]   - ūtrac)),   v = compute!(Field(trac[:v][i]   - v̄trac)))
    centervel = (u = compute!(Field(center[:u][i] - ūcenter)), v = compute!(Field(center[:v][i] - v̄center)))

    wenovort   = compute!(Field(VerticalVorticityField(weno, i)   - ζ̄weno))  
    tracvort   = compute!(Field(VerticalVorticityField(trac, i)   - ζ̄trac))  
    centervort = compute!(Field(VerticalVorticityField(center, i) - ζ̄center))

    set!(Eflucweno[i],   WenoNeverworld.Diagnostics.KineticEnergyField(wenovel))
    set!(Efluctrac[i],   WenoNeverworld.Diagnostics.KineticEnergyField(tracvel))
    set!(Efluccenter[i], WenoNeverworld.Diagnostics.KineticEnergyField(centervel))

    set!(Ωflucweno[i],   compute!(Field(wenovort^2)))
    set!(Ωfluctrac[i],   compute!(Field(tracvort^2)))
    set!(Ωfluccenter[i], compute!(Field(centervort^2)))
end

Eflucmeanweno1   = time_average(Eflucweno)
Eflucmeantrac1   = time_average(Efluctrac)
Eflucmeancenter1 = time_average(Efluccenter)

Ωflucmeanweno1   = time_average(Ωflucweno)
Ωflucmeantrac1   = time_average(Ωfluctrac)
Ωflucmeancenter1 = time_average(Ωfluccenter)

Eflucmeanweno2   = compute!(Field(Etotweno   - Ēweno))
Eflucmeantrac2   = compute!(Field(Etottrac   - Ētrac))
Eflucmeancenter2 = compute!(Field(Etotcenter - Ēcenter))

Ωflucmeanweno2   = compute!(Field(Ωtotweno   - Ω̄weno))
Ωflucmeantrac2   = compute!(Field(Ωtottrac   - Ω̄trac))
Ωflucmeancenter2 = compute!(Field(Ωtotcenter - Ω̄center))

IIĒweno   = compute!(Field(Integral(Ēweno,   dims = (1, 3))))
IIĒtrac   = compute!(Field(Integral(Ētrac,   dims = (1, 3))))
IIĒcenter = compute!(Field(Integral(Ēcenter, dims = (1, 3))))

IIEflucweno1   = compute!(Field(Integral(Eflucmeanweno1,   dims = (1, 3))))
IIEfluctrac1   = compute!(Field(Integral(Eflucmeantrac1,   dims = (1, 3))))
IIEfluccenter1 = compute!(Field(Integral(Eflucmeancenter1, dims = (1, 3))))

IIEflucweno2   = compute!(Field(Integral(Eflucmeanweno2,   dims = (1, 3))))
IIEfluctrac2   = compute!(Field(Integral(Eflucmeantrac2,   dims = (1, 3))))
IIEfluccenter2 = compute!(Field(Integral(Eflucmeancenter2, dims = (1, 3))))

IIEtotweno   = compute!(Field(Integral(Etotweno,   dims = (1, 3))))
IIEtottrac   = compute!(Field(Integral(Etottrac,   dims = (1, 3))))
IIEtotcenter = compute!(Field(Integral(Etotcenter, dims = (1, 3))))

IIΩ̄weno        = compute!(Field(Integral(Ω̄weno,   dims = (1, 3))))
IIΩ̄trac        = compute!(Field(Integral(Ω̄trac,   dims = (1, 3))))
IIΩ̄center      = compute!(Field(Integral(Ω̄center, dims = (1, 3))))

IIΩflucweno1   = compute!(Field(Integral(Ωflucmeanweno1,   dims = (1, 3))))
IIΩfluctrac1   = compute!(Field(Integral(Ωflucmeantrac1,   dims = (1, 3))))
IIΩfluccenter1 = compute!(Field(Integral(Ωflucmeancenter1, dims = (1, 3))))

IIΩflucweno2   = compute!(Field(Integral(Ωflucmeanweno2,   dims = (1, 3))))
IIΩfluctrac2   = compute!(Field(Integral(Ωflucmeantrac2,   dims = (1, 3))))
IIΩfluccenter2 = compute!(Field(Integral(Ωflucmeancenter2, dims = (1, 3))))

IIΩtotweno     = compute!(Field(Integral(Ωtotweno,   dims = (1, 3))))
IIΩtottrac     = compute!(Field(Integral(Ωtottrac,   dims = (1, 3))))
IIΩtotcenter   = compute!(Field(Integral(Ωtotcenter, dims = (1, 3))))

ϕ  = ynodes(IIEtotweno)
ϕf = ynodes(IIΩtotweno)

fig = Figure()
ax = Axis(fig[1, 1])

lines!(ax, ϕ, interior(IIEflucweno1  , 1, :, 1), color = :blue)
lines!(ax, ϕ, interior(IIEfluctrac1  , 1, :, 1), color = :green)
lines!(ax, ϕ, interior(IIEfluccenter1, 1, :, 1), color = :red)
lines!(ax, ϕ, interior(IIEflucweno2  , 1, :, 1), color = :blue , linestyle = :dash)
lines!(ax, ϕ, interior(IIEfluctrac2  , 1, :, 1), color = :green, linestyle = :dash)
lines!(ax, ϕ, interior(IIEfluccenter2, 1, :, 1), color = :red  , linestyle = :dash)

ax = Axis(fig[1, 2])

lines!(ax, ϕ, interior(IIĒweno  , 1, :, 1), color = :blue)
lines!(ax, ϕ, interior(IIĒtrac  , 1, :, 1), color = :green)
lines!(ax, ϕ, interior(IIĒcenter, 1, :, 1), color = :red)

ax = Axis(fig[1, 3])

lines!(ax, ϕ, interior(IIEtotweno  , 1, :, 1), color = :blue)
lines!(ax, ϕ, interior(IIEtottrac  , 1, :, 1), color = :green)
lines!(ax, ϕ, interior(IIEtotcenter, 1, :, 1), color = :red)

ax = Axis(fig[2, 1])

lines!(ax, ϕf, interior(IIΩflucweno1  , 1, :, 1), color = :blue)
lines!(ax, ϕf, interior(IIΩfluctrac1  , 1, :, 1), color = :green)
lines!(ax, ϕf, interior(IIΩfluccenter1, 1, :, 1), color = :red)
lines!(ax, ϕf, interior(IIΩflucweno2  , 1, :, 1), color = :blue , linestyle = :dash)
lines!(ax, ϕf, interior(IIΩfluctrac2  , 1, :, 1), color = :green, linestyle = :dash)
lines!(ax, ϕf, interior(IIΩfluccenter2, 1, :, 1), color = :red  , linestyle = :dash)

ax = Axis(fig[2, 2])

lines!(ax, ϕf, interior(IIΩ̄weno  , 1, :, 1), color = :blue)
lines!(ax, ϕf, interior(IIΩ̄trac  , 1, :, 1), color = :green)
lines!(ax, ϕf, interior(IIΩ̄center, 1, :, 1), color = :red)

ax = Axis(fig[2, 3])

lines!(ax, ϕf, interior(IIΩtotweno  , 1, :, 1), color = :blue)
lines!(ax, ϕf, interior(IIΩtottrac  , 1, :, 1), color = :green)
lines!(ax, ϕf, interior(IIΩtotcenter, 1, :, 1), color = :red)

display(fig)

IEtotweno   = compute!(Field(Integral(weno[:E][31],   dims = 3)))
IEtottrac   = compute!(Field(Integral(trac[:E][31],   dims = 3)))
IEtotcenter = compute!(Field(Integral(center[:E][31], dims = 3)))

fig = Figure()
ax  = Axis(fig[1, 1], title = "Depth integreted kinetic energy, WENO")
heatmap!(ax, log.(interior(IEtotweno,   :, :, 1) .+ 1e-20), colorrange = (0, log(1e3)), colormap = :magma)

ax  = Axis(fig[1, 2], title = "Depth integreted kinetic energy, rotated WENO")
heatmap!(ax, log.(interior(IEtottrac,   :, :, 1) .+ 1e-20), colorrange = (0, log(1e3)), colormap = :magma)

ax  = Axis(fig[1, 3], title = "Depth integreted kinetic energy, Centered")
heatmap!(ax, log.(interior(IEtotcenter, :, :, 1) .+ 1e-20), colorrange = (0, log(1e3)), colormap = :magma)

display(fig)

WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(center)
WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(weno)
WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(trac)

"""
Calculate energy (Ê) and enstrophy (Ω̂) spectra on a transect in the channel
"""

using WenoNeverworld.Diagnostics: average_spectra, hann_window

Êcenter = average_spectra(Efluccenter, Colon(), 80:80; k = 68)
Êweno   = average_spectra(Eflucweno,   Colon(), 80:80; k = 68)
Êtrac   = average_spectra(Efluctrac,   Colon(), 80:80; k = 68)

Ω̂center = average_spectra(Efluccenter, Colon(), 215:215; k = 68, windowing = hann_window)
Ω̂weno   = average_spectra(Eflucweno,   Colon(), 215:215; k = 68, windowing = hann_window)
Ω̂trac   = average_spectra(Efluctrac,   Colon(), 215:215; k = 68, windowing = hann_window)

NxE = length(xnodes(weno[:E][1]))
Nxζ = length(xnodes(weno[:ζ][1]))

fx = Êcenter.freq[2:end]
fxE = fftfreq(NxE, 1 / Δxᶜᶜᶜ(1, 1, 1, weno[:E].grid))[1:Int(NxE ÷ 2)] .* 1e3 # in 1/km
fxζ = fftfreq(Nxζ, 1 / Δxᶠᶠᶜ(1, 1, 1, weno[:ζ].grid))[1:Int(Nxζ ÷ 2)] .* 1e3 # in 1/km

fig = Figure()
ax = Axis(fig[1, 1], title = "energy spectra", yscale = log10, xscale = log10)
lines!(ax, fxE[2:end], Êcenter.spec[2:end], color = :red,   label = "center momentum scheme")
lines!(ax, fxE[2:end], Êweno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
lines!(ax, fxE[2:end], Êtrac.spec[2:end],   color = :green,  label = "WENO momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end] ./ Êweno.spec[2] .* Êcenter.spec[2], color = :green, linestyle = :dashdot)
lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

ax = Axis(fig[1, 2], title = "enstrophy spectra", yscale = log10, xscale = log10)
lines!(ax, fxE[2:end], Ω̂center.spec[2:end], color = :red,   label = "center momentum scheme")
lines!(ax, fxE[2:end], Ω̂weno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
lines!(ax, fxE[2:end], Ω̂trac.spec[2:end],   color = :green, label = "WENO momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end] ./ Êweno.spec[2] .* Êcenter.spec[2], color = :green, linestyle = :dashdot)
lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

display(fig)
