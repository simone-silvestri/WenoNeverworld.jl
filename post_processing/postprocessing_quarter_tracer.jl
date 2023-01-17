using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Statistics: mean
using GLMakie

"""
Time-averaged variables
"""

# variables = ("b", ) # "u", "v")

# centerb_avg = all_fieldtimeseries("../UpwindVectorInvariantSchemes.jl/quarter_centered/global_upwinding_averages.jld2"; variables)
# leithb_avg  = all_fieldtimeseries("../UpwindVectorInvariantSchemes.jl/quarter_leith/global_upwinding_averages.jld2"; variables)
# wenob_avg   = all_fieldtimeseries("../UpwindVectorInvariantSchemes.jl/quarter/global_upwinding_averages.jld2"; variables)
# eightb_avg  = all_fieldtimeseries("../UpwindVectorInvariantSchemes.jl/eighth/global_upwinding_eight_averages.jld2"; variables)

# # centerb_avg = all_fieldtimeseries("centered_check", "../UpwindVectorInvariantSchemes.jl/tracer_experiment/"; variables, checkpointer = true)
# # leithb_avg  = all_fieldtimeseries("centered_leith", "../UpwindVectorInvariantSchemes.jl/tracer_experiment/"; variables, checkpointer = true)
# # wenob_avg   = all_fieldtimeseries("weno_quarter",   "../UpwindVectorInvariantSchemes.jl/tracer_experiment/"; variables, checkpointer = true)
# # eightb_avg  = all_fieldtimeseries("weno_eight",     "../UpwindVectorInvariantSchemes.jl/tracer_experiment/"; variables, checkpointer = true)

# times  = centerb_avg[:b].times[end-33:end]
# times8 =  eightb_avg[:b].times[end-33:end]

# centerb_avg = limit_timeseries!(centerb_avg, times)
# leithb_avg  = limit_timeseries!(leithb_avg,  times)
# wenob_avg   = limit_timeseries!(wenob_avg,   times)
# eightb_avg  = limit_timeseries!(eightb_avg,  times8)

# # # @info "computing integrated heat content..."
# # # centerb_heat = WenoNeverworld.Diagnostics.heat_content(centerb_avg[:b])
# # # wenob_heat   = WenoNeverworld.Diagnostics.heat_content(wenob_avg[:b])

# @info "computing time averaged buoyancy fields..."
# bwmean = time_average(wenob_avg[:b]) 
# bcmean = time_average(centerb_avg[:b])
# blmean = time_average(leithb_avg[:b])
# b8mean = time_average(eightb_avg[:b])

# @info "computing zonally and time averaged buoyancy fields..."
# Bw = compute!(Field(Average(bwmean, dims = 1)))
# Bc = compute!(Field(Average(bcmean, dims = 1)))
# Bl = compute!(Field(Average(blmean, dims = 1)))
# B8 = compute!(Field(Average(b8mean, dims = 1)))

# @info "plotting"
# ϕ = ynodes(wenob_avg[:b][1])
# z = znodes(wenob_avg[:b][1])

# ϕc = ynodes(centerb_avg[:b][1])
# zc = znodes(centerb_avg[:b][1])

# ϕ8 = ynodes(eightb_avg[:b][1])
# z8 = znodes(eightb_avg[:b][1])

# fig = Figure(resolution = (2400 * 5, 500 * 5))
# ax  = Axis(fig[1, 1], title = "C4 - W8 comparison", xlabel = L"\phi", ylabel = L"z",
#            xgridvisible = false, ygridvisible = false)
# cr = contourf!(ax, ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), colormap = :thermometer)
# contour!(ax, ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :dash, linewidth = 5, color = :black)
# contour!(ax, ϕ, z[6:end] ./ 1000,  interior(Bc, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :solid, linewidth = 10, color = :black)
# ylims!(ax, (z[6], z[end])./ 1000)
# xlims!(ax, (-70, 0))

# ax  = Axis(fig[1, 2], title = "L4 - W8 comparison", xlabel = L"\phi",
#            xgridvisible = false, ygridvisible = false)
# cr = contourf!(ax,      ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), colormap = :thermometer)
# contour!(ax, ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :dash, linewidth = 5, color = :black)
# contour!(ax, ϕ, z[6:end] ./ 1000,  interior(Bl, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :solid, linewidth = 10, color = :black)
# ylims!(ax, (z[6], z[end])./ 1000)
# xlims!(ax, (-70, 0))

# ax  = Axis(fig[1, 3], title = "W4 - W8 comparison", xlabel = L"\phi",
#            xgridvisible = false, ygridvisible = false)
# cr = contourf!(ax,      ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), colormap = :thermometer)
# contour!(ax, ϕ8, z8[6:end] ./ 1000, interior(B8, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :dash, linewidth = 5, color = :black)
# contour!(ax, ϕ, z[6:end] ./ 1000,  interior(Bw, 1, :, 6:69), levels = (0.0000000001:0.002:0.06), linestyle = :solid,  linewidth = 10, color = :black)
# ylims!(ax, (z[6], z[end])./ 1000)
# xlims!(ax, (-70, 0))

# cb = Colorbar(fig[1, 4], cr)

# save("buoyancy.png", fig)
# using CairoMakie
# CairoMakie.save("buoyancy.eps", fig)

# using WenoNeverworld.Diagnostics: DensityField

# @info "computing equilibrium height z★..."

# function compute_rpe_density(var)
#     ze = WenoNeverworld.Diagnostics.calculate_z★_diagnostics(var[:b])

#     εe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)
#     αe = FieldTimeSeries{Center, Center, Center}(ze.grid, ze.times)

#     zfield = WenoNeverworld.Diagnostics.HeightField(ze.grid)

#     @info "computing resting and available potential energy density..."
#     for t in 1:length(ze.times)
#         @info "doing time $t"
#         ρ = DensityField(var[:b][t])
#         set!(εe[t], compute!(Field(ze[t] * ρ)))
#         set!(αe[t], compute!(Field((- zfield - ze[t]) * ρ)))
#     end

#     εeavg = mean(compute!(Field((εe[length(ze.times)] - εe[1]))), dims = 1)

#     return (; ze, εe, αe, εeavg)
# end

# stc = compute_rpe_density(centerb_avg)
# stl = compute_rpe_density(leithb_avg )
# stw = compute_rpe_density(wenob_avg  )
# ste = compute_rpe_density(eightb_avg )

# function calculate_RPE(st)
#     RPE = Float64[]
#     vol = VolumeField(st.ze.grid)

#     for t in 1:length(st.ze.times)
#         @info "doing time $t"
#         push!(RPE, sum(compute!(Field(st.εe[t] * vol))))
#     end

#     return RPE
# end

# function calculate_APE(st)
#     APE = Float64[]
#     vol = VolumeField(st.ze.grid)

#     for t in 1:length(st.ze.times)
#         @info "doing time $t"
#         push!(APE, sum(compute!(Field(st.αe[t] * vol))))
#     end

#     return APE
# end

# function calculate_KE(var)
#     KE  = Float64[]
#     vol = VolumeField(var[:u].grid)

#     @info "computing resting and available potential energy density..."
#     for t in 1:length(var[:u].times)
#         @info "doing time $t"
#         v = var[:v][t]
#         u = var[:u][t]
#         ke = compute!(Field(@at (Center, Center, Center) u^2 + v^2))

#         push!(KE, sum(compute!(Field(ke * vol))) * 0.5)
#     end

#     return KE
# end

# RPEc = calculate_RPE(stc)
# RPEl = calculate_RPE(stl)
# RPEw = calculate_RPE(stw)
# RPEe = calculate_RPE(ste)

# APEc = calculate_APE(stc)
# APEl = calculate_APE(stl)
# APEw = calculate_APE(stw)
# APEe = calculate_APE(ste)

# KEc = calculate_KE(centerb_avg)
# KEl = calculate_KE( leithb_avg)
# KEw = calculate_KE(  wenob_avg)
# KEe = calculate_KE( eightb_avg)

# fig = Figure()
# ax  = Axis(fig[1, 1], title = "Zonally averaged buoyancy")
# contourf!(ax, ϕ, z[6:end], interior(εcavg, 1, :, 6:69); levels)
# contour!(ax,  ϕ, z[6:end], interior(εwavg, 1, :, 6:69); levels, color = :black)

# times4 =  wenob_avg[:u].times .-  wenob_avg[:u].times[1]
# times8 = eightb_avg[:u].times .- eightb_avg[:u].times[1]

# color1 = :deepskyblue
# color2 = :orange1
# color3 = :firebrick2

# fig = Figure()
# ax  = Axis(fig[1, 1], xgridvisible = false, ygridvisible = false, ylabel = L"(RPE - RPE(0)) / RPE(0)", 
#            xlabel = "months", title = L"RPE = \int_V \rho z^\star dV",
#            xticks = ([0, 4, 8, 12], [L"0", L"4", L"8", L"12"]),
#            yticks = ([0, 5e-5, 1e-4, 1.5e-4], [L"0.0", L"5\cdot 10^{-5}", L"1\cdot 10^{-4}", L"1.5\cdot 10^{-4}"]))
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (RPEe .- RPEe[1]) ./ RPEe[1], color = :black, linewidth = 2, label = L"W8", linestyle = :dash)
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (RPEl .- RPEl[1]) ./ RPEl[1], color = color1, linewidth = 2, label = L"L4")
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (RPEc .- RPEc[1]) ./ RPEc[1], color = color2, linewidth = 2, label = L"C4")
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (RPEw .- RPEw[1]) ./ RPEw[1], color = color3, linewidth = 2, label = L"W4")
# xlims!(ax, (0, 12))
# ylims!(ax, (0.0, 1.75e-4))
# axislegend(ax, position = :lt)

# ax  = Axis(fig[1, 2],  xgridvisible = false, ygridvisible = false, ylabel = L"(APE - APE(0)) / APE(0)",
#            xlabel = "months", title = L"APE = \int_V \rho (z - z^\star) dV",
#            xticks = ([0, 4, 8, 12], [L"0", L"4", L"8", L"12"]),
#            yticks = ([-1.5e-3, -1e-3, -5e-4, 0.0], [L"-1.5\cdot 10^{-3}", L"-1\cdot 10^{-3}", L"-5\cdot 10^{-4}", L"0.0"]))
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (APEe .- APEe[1]) ./ APEe[1], color = :black, linewidth = 2, label = L"W8", linestyle = :dash)
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (APEl .- APEl[1]) ./ APEl[1], color = color1, linewidth = 2, label = L"L4")
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (APEc .- APEc[1]) ./ APEc[1], color = color2, linewidth = 2, label = L"C4")
# lines!(ax, times4 ./ 60 ./ 24 ./ 30, (APEw .- APEw[1]) ./ APEw[1], color = color3, linewidth = 2, label = L"W4")
# xlims!(ax, (0, 12))
# ylims!(ax, (-1.75e-3, 0.0))

# ax  = Axis(fig[1, 3], title = "(KE - KE(0)")
# lines!(ax, times4, (KEe .- KEe[1]), color = :black, linewidth = 2)
# lines!(ax, times4, (KEl .- KEl[1]), color = color1, linewidth = 2)
# lines!(ax, times4, (KEc .- KEc[1]), color = color2, linewidth = 2)
# lines!(ax, times4, (KEw .- KEw[1]), color = color3, linewidth = 2)
