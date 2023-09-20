using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using JLD2
using Statistics
using GLMakie

output_dir = "/storage4/WenoNeverworldData/"
output_prefix = output_dir * "weno_eight_surface.jld2"
jlfile = jldopen(output_prefix)
utimeseries = jlfile["timeseries"]["u"]
vtimeseries = jlfile["timeseries"]["v"]
btimeseries = jlfile["timeseries"]["b"]
ukeys = keys(utimeseries)[2:end]
fig = Figure(resolution=(600, 1200))
ax11 = Axis(fig[1, 1])
ax12 = Axis(fig[2, 1])
ax13 = Axis(fig[3, 1])
ufield = utimeseries[ukeys[end]][:, :, 1]
vfield = vtimeseries[ukeys[end]][:, :, 1]
bfield = btimeseries[ukeys[end]][:, :, 1]
heatmap!(ax11, ufield, colormap=:balance, colorrange=(-3, 3))
heatmap!(ax12, vfield, colormap=:balance, colorrange=(-3, 3))
heatmap!(ax13, bfield, colormap=:plasma, colorrange=extrema(bfield))
display(fig)
close(jlfile)

##
output_dir = "/storage4/WenoNeverworldData/"
output_prefix = output_dir * "weno_eight_snapshots.jld2"
jlfile = jldopen(output_prefix)
utimeseries = jlfile["timeseries"]["u"]
vtimeseries = jlfile["timeseries"]["v"]
btimeseries = jlfile["timeseries"]["b"]
imgrid = jlfile["serialized"]["grid"]
z = imgrid.underlying_grid.zᵃᵃᶜ[1:69]
λᶜᵃᵃ = imgrid.underlying_grid.λᶜᵃᵃ[1:512]
φᵃᶠᵃ = imgrid.underlying_grid.φᵃᶠᵃ[1:401]
φᵃᶜᵃ = imgrid.underlying_grid.φᵃᶜᵃ[1:400]

ukeys = keys(utimeseries)[2:end]
fig2 = Figure(resolution=(2000, 1000))
ax11 = Axis(fig2[1, 1])
ax12 = Axis(fig2[2, 1])
ax13 = Axis(fig2[3, 1])
# ufield = mean(utimeseries[ukeys[end]], dims=1)[1, :, :]
# vfield = mean(vtimeseries[ukeys[end]], dims=1)[1, :, :]
# bfield = mean(btimeseries[ukeys[end]], dims=1)[1, :, :]
sl_x = Slider(fig2[4, 1], range=1:512, startvalue=3)
obs = sl_x.value
ufield = @lift utimeseries[ukeys[end]][$obs, :, :]
vfield = @lift vtimeseries[ukeys[end]][$obs, :, :]
bfield = @lift btimeseries[ukeys[end]][$obs, :, :]
heatmap!(ax11, φᵃᶜᵃ, z, ufield, colormap=:balance, colorrange=(-3, 3))
heatmap!(ax12, φᵃᶠᵃ, z, vfield, colormap=:balance, colorrange=(-3, 3))
heatmap!(ax13, φᵃᶜᵃ, z, bfield, colormap=:plasma, colorrange=(0, 0.05))
contour!(ax13, φᵃᶜᵃ, z, bfield, levels=20, linewidth=5, color=:black)
display(fig2)
close(jlfile)

##
output_dir = "/storage2/WenoNeverworldData/"
output_prefix = output_dir * "weno_two_snapshots.jld2"
jlfile_2 = jldopen(output_prefix)
utimeseries_2 = jlfile_2["timeseries"]["u"]
vtimeseries_2 = jlfile_2["timeseries"]["v"]
btimeseries_2 = jlfile_2["timeseries"]["b"]
tkeys_2 = keys(utimeseries_2)[2:end]
imgrid = jlfile_2["serialized"]["grid"]
λᶜᵃᵃ_2_size, φᵃᶜᵃ_2_size, z_size_2 = size(btimeseries_2[tkeys_2[end]])
z_2 = imgrid.underlying_grid.zᵃᵃᶜ[1:z_size_2]
λᶜᵃᵃ_2 = imgrid.underlying_grid.λᶜᵃᵃ[1:λᶜᵃᵃ_2_size]
φᵃᶠᵃ_2 = imgrid.underlying_grid.φᵃᶠᵃ[1:φᵃᶜᵃ_2_size+1]
φᵃᶜᵃ_2 = imgrid.underlying_grid.φᵃᶜᵃ[1:φᵃᶜᵃ_2_size]

output_prefix = output_dir * "weno_four_snapshots.jld2"
jlfile_4 = jldopen(output_prefix)
utimeseries_4 = jlfile_4["timeseries"]["u"]
vtimeseries_4 = jlfile_4["timeseries"]["v"]
btimeseries_4 = jlfile_4["timeseries"]["b"]
tkeys_4 = keys(utimeseries_4)[2:end]
imgrid = jlfile_4["serialized"]["grid"]
λᶜᵃᵃ_4_size, φᵃᶜᵃ_4_size, z_size_4 = size(btimeseries_4[tkeys_4[end]])
z_4 = imgrid.underlying_grid.zᵃᵃᶜ[1:z_size_4]
λᶜᵃᵃ_4 = imgrid.underlying_grid.λᶜᵃᵃ[1:λᶜᵃᵃ_4_size]
φᵃᶠᵃ_4 = imgrid.underlying_grid.φᵃᶠᵃ[1:φᵃᶜᵃ_4_size+1]
φᵃᶜᵃ_4 = imgrid.underlying_grid.φᵃᶜᵃ[1:φᵃᶜᵃ_4_size]

output_prefix = output_dir * "weno_sixteenth_snapshots.jld2"
jlfile_16 = jldopen(output_prefix)
utimeseries_16 = jlfile_16["timeseries"]["u"]
vtimeseries_16 = jlfile_16["timeseries"]["v"]
btimeseries_16 = jlfile_16["timeseries"]["b"]
tkeys_16 = keys(utimeseries_16)[2:end]
imgrid = jlfile_16["serialized"]["grid"]
λᶜᵃᵃ_16_size, φᵃᶜᵃ_16_size, z_size_16 = size(btimeseries_16[tkeys_16[end]])
z_16 = imgrid.underlying_grid.zᵃᵃᶜ[1:z_size_16]
λᶜᵃᵃ_16 = imgrid.underlying_grid.λᶜᵃᵃ[1:λᶜᵃᵃ_16_size]
φᵃᶠᵃ_16 = imgrid.underlying_grid.φᵃᶠᵃ[1:φᵃᶜᵃ_16_size+1]
φᵃᶜᵃ_16 = imgrid.underlying_grid.φᵃᶜᵃ[1:φᵃᶜᵃ_16_size]

ukeys = keys(utimeseries)[2:end]
fig3 = Figure(resolution=(2000, 1000))
ax11 = Axis(fig3[1, 1])
ax12 = Axis(fig3[2, 1])
ax13 = Axis(fig3[1, 2])
ax14 = Axis(fig3[2, 2])
sl_x = Slider(fig3[3, 1:2], range=1:128, startvalue=3)
obs = sl_x.value
bfield_2 = @lift btimeseries_2[tkeys_2[end]][$obs, :, :]
bfield_4 = @lift btimeseries_4[tkeys_4[end]][$obs*2, :, :]
bfield   = @lift btimeseries[ukeys[end]][$obs*4, :, :]
bfield_16   = @lift btimeseries_16[tkeys_16[end]][$obs*8, :, :]
endval = 1.0
bcontours = (exp.(range(0, endval, length=20)) .- 1) / (exp(endval) - 1) * 0.05 
heatmap!(ax11, φᵃᶜᵃ_2, z_2, bfield_2, colormap=:plasma, colorrange=(0, 0.05))
contour!(ax11, φᵃᶜᵃ_2, z_2, bfield_2, levels=bcontours, linewidth=5, color=:black)
heatmap!(ax12, φᵃᶜᵃ_4, z_4, bfield_4, colormap=:plasma, colorrange=(0, 0.05))
contour!(ax12, φᵃᶜᵃ_4, z_4, bfield_4, levels=bcontours, linewidth=5, color=:black)
heatmap!(ax13, φᵃᶜᵃ, z, bfield, colormap=:plasma, colorrange=(0, 0.05))
contour!(ax13, φᵃᶜᵃ, z, bfield, levels=bcontours, linewidth=5, color=:black)
heatmap!(ax14, φᵃᶜᵃ_16, z_16, bfield_16, colormap=:plasma, colorrange=(0, 0.05))
contour!(ax14, φᵃᶜᵃ_16, z_16, bfield_16, levels=bcontours, linewidth=5, color=:black)
display(fig3)

#=
close(jlfile_2)
close(jlfile_4)
close(jlfile_16)
=#