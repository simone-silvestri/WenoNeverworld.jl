using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Operators
using GLMakie

"""
Time-averaged variables
"""

variables = ("b", )

center_avg = all_fieldtimeseries("files_four_centered_new_bathy/neverworld_quarter_centered_averages.jld2"; variables)
center_avg = limit_timeseries!(center_avg, center_avg[:b].times[end-50:end])

weno_avg = all_fieldtimeseries("files_four_new_bathy/neverworld_quarter_averages.jld2"; variables)
weno_avg = limit_timeseries!(weno_avg, weno_avg[:b].times[end-50:end])

trac_avg = all_fieldtimeseries("files_four_isopycnal_new_bathy/neverworld_quarter_isopycnal_averages.jld2"; variables)
trac_avg = limit_timeseries!(trac_avg, trac_avg[:b].times[end-50:end])

@info "computing integrated heat content..."
center_heat = WenoNeverworld.Diagnostics.heat_content(center_avg2[:b])
weno_heat   = WenoNeverworld.Diagnostics.heat_content(weno_avg2[:b])
trac_heat   = WenoNeverworld.Diagnostics.heat_content(trac_avg2[:b])

@info "computing time averaged buoyancy fields..."
bwmean = time_average(weno_avg[:b]) 
bcmean = time_average(center_avg[:b])
btmean = time_average(trac_avg[:b])  

@info "computing zonally and time averaged buoyancy fields..."
Bw = compute!(Field(Average(bwmean, dims = 1)))
Bc = compute!(Field(Average(bcmean, dims = 1)))
Bt = compute!(Field(Average(btmean, dims = 1)))

ϕ = ynodes(weno_avg[:b][1])
z = znodes(weno_avg[:b][1])

fig = Figure()
ax  = Axis(fig[1, 1], title = "Zonally averaged buoyancy")
contourf!(ax, ϕ, z, interior(Bw, 1, :, :), levels = (0:0.002:0.06))
contour!(ax,  ϕ, z, interior(Bc, 1, :, :), levels = (0:0.002:0.06), color = :black)
contour!(ax,  ϕ, z, interior(Bt, 1, :, :), levels = (0:0.002:0.06), color = :black, linestyle = :dash)

display(fig)

using WenoNeverworld.Diagnostics: DensityField

@info "computing equilibrium height z★..."
zc = WenoNeverworld.Diagnostics.calculate_z★_diagnostics(center_avg[:b])
zw = WenoNeverworld.Diagnostics.calculate_z★_diagnostics(weno_avg[:b])
zt = WenoNeverworld.Diagnostics.calculate_z★_diagnostics(trac_avg[:b])

εc = FieldTimeSeries{Center, Center, Center}(zc.grid, zc.times)
εw = FieldTimeSeries{Center, Center, Center}(zw.grid, zw.times)
εt = FieldTimeSeries{Center, Center, Center}(zt.grid, zt.times)

@info "computing resting potential energy density..."
for t in 1:length(zc.times)
    set!(εc[t], compute!(Field(zc[t] * DensityField(center_avg[:b][t]))))
    set!(εw[t], compute!(Field(zw[t] * DensityField(weno_avg[:b][t]))))
    set!(εt[t], compute!(Field(zt[t] * DensityField(trac_avg[:b][t]))))
end

vol = VolumeField(εc.grid)

RPEcini = sum(compute!(Field(εc[1] * vol)))
RPEwini = sum(compute!(Field(εw[1] * vol)))
RPEtini = sum(compute!(Field(εt[1] * vol)))

RPEc = zeros(length(εc.times))
RPEw = zeros(length(εw.times))
RPEt = zeros(length(εt.times))

@info "computing integral resting potential energy..."
for t in 2:length(zc.times) 
    RPEc[t] = sum(compute!(Field(εc[t] * vol))) - RPEcini
    RPEw[t] = sum(compute!(Field(εw[t] * vol))) - RPEwini
    RPEt[t] = sum(compute!(Field(εt[t] * vol))) - RPEtini
end