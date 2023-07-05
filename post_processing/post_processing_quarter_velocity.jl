using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Operators
using GLMakie
using FFTW 

using Statistics: mean

dir = "./limited_results/"

dict_or_file = :dict

if dict_or_file == :file
    """
    Instantaneous variables
    """

    variables = ("u", "v")

    center = all_fieldtimeseries(dir * "quarter_centered/global_upwinding_snapshots.jld2"; variables);
    weno   = all_fieldtimeseries(dir * "quarter/global_upwinding_snapshots.jld2"; variables);
    leith  = all_fieldtimeseries(dir * "quarter_leith/global_upwinding_snapshots.jld2"; variables);
    eight  = all_fieldtimeseries(dir * "eighth/global_upwinding_eight_snapshots.jld2"; variables);

    times4 = center[:u].times[end-20:end]
    times8 =  eight[:u].times[end-20:end]

    function calculate_noise(v)
        Nx, Ny, Nz = size(v)

        nx, ny, nz = (Nx - 100, Ny - 100, Nz - 20)

        χ = zeros(nx, ny, nz)

        for i in 51:Nx-50, j in 51:Ny-50, k in 21:Nz
            vm = 0.5 * v[i, j, k] + 0.25 * (v[i-1, j, k] + v[i+1, j, k])
            χ[i-50, j-50, k-20] = abs(v[i, j, k] - vm)
        end

        @info "done calculating noise"
        return χ
    end

    # finalc = length(center[:v].times)
    # finalw = length(  weno[:v].times)
    # finall = length( leith[:v].times)
    # finale = length( eight[:v].times)

    # χc = calculate_noise(center[:v][finalc])
    # χw = calculate_noise(  weno[:v][finalw])
    # χl = calculate_noise( leith[:v][finall])
    # χe = calculate_noise( eight[:v][finale])

    # χ̄c = dropdims(mean(χc, dims = (1, 3)), dims = (1, 3))
    # χ̄w = dropdims(mean(χw, dims = (1, 3)), dims = (1, 3))
    # χ̄l = dropdims(mean(χl, dims = (1, 3)), dims = (1, 3))
    # χ̄e = dropdims(mean(χe, dims = (1, 3)), dims = (1, 3))

    center = limit_timeseries!(center, times4) 
    weno   = limit_timeseries!(weno  , times4) 
    leith  = limit_timeseries!(leith , times4) 
    eight  = limit_timeseries!(eight , times8) 

    @info "finished loading instantaneous files"

    """
    Time-averaged variables
    """

    variables = ("u", "v", "u2", "v2", "ζ", "ζ2")

    center_avg = all_fieldtimeseries(dir * "quarter_centered/global_upwinding_averages.jld2"; variables);
    weno_avg   = all_fieldtimeseries(dir * "quarter/global_upwinding_averages.jld2"; variables);
    leith_avg  = all_fieldtimeseries(dir * "quarter_leith/global_upwinding_averages.jld2"; variables);
    eight_avg  = all_fieldtimeseries(dir * "eighth/global_upwinding_eight_averages.jld2"; variables);

    center_avg = limit_timeseries!(center_avg, times4) 
    weno_avg   = limit_timeseries!(weno_avg  , times4) 
    leith_avg  = limit_timeseries!(leith_avg , times4) 
    eight_avg  = limit_timeseries!(eight_avg , times8) 

    @info "Finished loading averaged files"

else
    """
    Instantaneous variables
    """

    center = jldopen(dir * "center.jld2")["vars"]
    weno   = jldopen(dir * "weno.jld2")["vars"]
    leith  = jldopen(dir * "leith.jld2")["vars"]
    eight  = jldopen(dir * "eight.jld2")["vars"]

    """
    Time-averaged variables
    """
    
    center_avg = jldopen(dir * "center_avg.jld2")["vars"]
    weno_avg   = jldopen(dir * "weno_avg.jld2")["vars"]
    leith_avg  = jldopen(dir * "leith_avg.jld2")["vars"]
    eight_avg  = jldopen(dir * "eight_avg.jld2")["vars"]
end

"""
Mean Kinetic energy and Enstrophy
"""

using WenoNeverworld.Diagnostics: VerticalVorticityField

WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(center)
WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(weno  )
WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(leith )
WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(eight )

function mean_kinetic_energy_and_enstrophy(var, var_avg)
    ū = time_average(var_avg[:u])
    v̄ = time_average(var_avg[:v])
    ζ̄ = time_average(var_avg[:ζ])

    WenoNeverworld.Diagnostics.add_kinetic_energy_from_timeseries!(var_avg)

    Etot = time_average(var_avg[:E])
    Ωtot = time_average(var_avg[:ζ2])

    Ē = WenoNeverworld.Diagnostics.KineticEnergyField((u = ū, v = v̄))
    Ω̄ = compute!(Field(ζ̄^2))

    Efluc = FieldTimeSeries{Center, Center, Center}(var[:u].grid, var[:u].times)
    Ωfluc = FieldTimeSeries{Face, Face, Center}(var[:u].grid, var[:u].times)

    for i in 1:length(var[:u].times)
        @info "doing time $i"
        vel  = (u = compute!(Field(var[:u][i] - ū)), v = compute!(Field(var[:v][i] - v̄)))
        vort = compute!(Field(VerticalVorticityField(var, i) - ζ̄))  
        set!(Efluc[i], WenoNeverworld.Diagnostics.KineticEnergyField(vel))
        set!(Ωfluc[i], compute!(Field(vort^2)))
    end

    Eflucmean1   = time_average(Efluc)
    Ωflucmean1   = time_average(Ωfluc)
    Eflucmean2   = compute!(Field(Etot - Ē))
    Ωflucmean2   = compute!(Field(Ωtot - Ω̄))

    IIĒ        = compute!(Field(Integral(Ē,          dims = (1, 3))))
    IIEfluc1   = compute!(Field(Integral(Eflucmean1, dims = (1, 3))))
    IIEfluc2   = compute!(Field(Integral(Eflucmean2, dims = (1, 3))))
    IIEtot     = compute!(Field(Integral(Etot,       dims = (1, 3))))
    IIΩ̄        = compute!(Field(Integral(Ω̄,          dims = (1, 3))))
    IIΩfluc1   = compute!(Field(Integral(Ωflucmean1, dims = (1, 3))))
    IIΩfluc2   = compute!(Field(Integral(Ωflucmean2, dims = (1, 3))))
    IIΩtot     = compute!(Field(Integral(Ωtot,       dims = (1, 3))))

    return (; Etot, Ωtot, Ē, Ω̄, Eflucmean1, Ωflucmean1, IIEtot, IIĒ, IIEfluc1, IIEfluc2, IIΩtot, IIΩ̄, IIΩfluc1, IIΩfluc2)
end

function barotropic_and_baroclinic_energy(var)

    Ebt  = FieldTimeSeries{Center, Center, Center}(var[:u].grid, var[:u].times)
    Ebtf = FieldTimeSeries{Center, Center, Center}(var[:u].grid, var[:u].times)

    for i in 1:length(var[:u].times)
        @info "doing time $i"
        ubt  = compute!(Field(Integral(var[:u][i], dims = 3)))
        vbt  = compute!(Field(Integral(var[:v][i], dims = 3)))
        vel  = (u = ubt, v = vbt)
        Eint = compute!(Field(Integral(var[:E][i], dims = 3)))
        set!(Ebt[i],  WenoNeverworld.Diagnostics.KineticEnergyField(vel))
        set!(Ebtf[i], compute!(Field(Ebt[i] / Eint)))
    end

    return compute!(Field(Average(time_average(Ebtf), dims = 1)))
end

btc = barotropic_and_baroclinic_energy(center)
btw = barotropic_and_baroclinic_energy(weno)
btl = barotropic_and_baroclinic_energy(leith)
bte = barotropic_and_baroclinic_energy(eight)

center_avg = mean_kinetic_energy_and_enstrophy(center, center_avg)
weno_avg   = mean_kinetic_energy_and_enstrophy(weno  , weno_avg  )
leith_avg  = mean_kinetic_energy_and_enstrophy(leith , leith_avg )
eight_avg  = mean_kinetic_energy_and_enstrophy(eight , eight_avg )

ϕ   = ynodes(leith[:E][1])
ϕf  = ynodes(leith[:ζ][1])
ϕ8  = ynodes(eight[:E][1])
ϕ8f = ynodes(eight[:ζ][1])

fig = Figure()
ax = Axis(fig[1, 1])

lines!(ax, ϕ8, interior(eight_avg.IIEfluc1,  1, :, 1), color = :black, linestyle = :dash)
lines!(ax, ϕ,  interior(center_avg.IIEfluc1, 1, :, 1), color = :blue,  linestyle = :dash)
lines!(ax, ϕ,  interior(weno_avg.IIEfluc1,   1, :, 1), color = :red,   linestyle = :dash)
lines!(ax, ϕ,  interior(leith_avg.IIEfluc1,  1, :, 1), color = :green, linestyle = :dash)

lines!(ax, ϕ8, interior(eight_avg.IIEfluc1,  1, :, 1) .+ interior(eight_avg.IIĒ,  1, :, 1), color = :black)
lines!(ax, ϕ,  interior(center_avg.IIEfluc1, 1, :, 1) .+ interior(center_avg.IIĒ, 1, :, 1), color = :blue )
lines!(ax, ϕ,  interior(weno_avg.IIEfluc1,   1, :, 1) .+ interior(weno_avg.IIĒ,   1, :, 1), color = :red  )
lines!(ax, ϕ,  interior(leith_avg.IIEfluc1,  1, :, 1) .+ interior(leith_avg.IIĒ,  1, :, 1), color = :green)

ax = Axis(fig[1, 2])

lines!(ax, ϕ8f, interior(eight_avg.IIΩfluc1,  1, :, 1), color = :black, linestyle = :dash)
lines!(ax, ϕf,  interior(center_avg.IIΩfluc1, 1, :, 1), color = :blue,  linestyle = :dash)
lines!(ax, ϕf,  interior(weno_avg.IIΩfluc1,   1, :, 1), color = :red,   linestyle = :dash)
lines!(ax, ϕf,  interior(leith_avg.IIΩfluc1,  1, :, 1), color = :green, linestyle = :dash)

lines!(ax, ϕ8f, interior(eight_avg.IIΩtot,  1, :, 1), color = :black)
lines!(ax, ϕf,  interior(center_avg.IIΩtot, 1, :, 1), color = :blue)
lines!(ax, ϕf,  interior(weno_avg.IIΩtot,   1, :, 1), color = :red)
lines!(ax, ϕf,  interior(leith_avg.IIΩtot,  1, :, 1), color = :green)

display(fig)

# CairoMakie.save("tot_energy.eps", fig)

# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(center)
# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(weno)
# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(leith)
# WenoNeverworld.Diagnostics.add_kinetic_energy_and_vorticity_timeseries!(eight)

# λ  = xnodes(leith[:E][1])
# λ8 = xnodes(eight[:E][1])

# IEtotweno   = compute!(Field(Integral(  weno[:E][length(times4)], dims = 3)))
# IEtotcenter = compute!(Field(Integral(center[:E][length(times4)], dims = 3)))
# IEtotleith  = compute!(Field(Integral( leith[:E][length(times4)], dims = 3)))
# IEtoteight  = compute!(Field(Integral( eight[:E][length(times8)], dims = 3)))


# fig = Figure(resolution = (7000, 2000))
# ax  = Axis(fig[1, 1], title = "depth-integrated TKE, weno scheme", xlabel = L"\lambda", ylabel = L"\phi")
# # heatmap!(ax, λ, ϕ, log10.(interior(IEtotweno,   :, :, 1) .+ 1e-20), colorrange = (0, log10(1e3)), colormap = :magma)
# heatmap!(ax, λ, ϕ, log10.(interior(weno[:E][1], :, :, 65) .+ 1e-20), colorrange = (-5, 0), colormap = :magma, interpolate = true)

# ax = Axis(fig[1, 2], title = "depth-integrated TKE, centered scheme", xlabel = L"\lambda", ylabel = L"\phi")
# heatmap!(ax, λ, ϕ, log10.(interior(center[:E][1], :, :, 65) .+ 1e-20), colorrange = (-5, 0), colormap = :magma, interpolate = true)

# ax = Axis(fig[1, 3], title = "depth-integrated TKE, leith scheme", xlabel = L"\lambda", ylabel = L"\phi")
# heatmap!(ax, λ, ϕ, log10.(interior(leith[:E][1], :, :, 65) .+ 1e-20), colorrange = (-5, 0), colormap = :magma, interpolate = true)

# ax = Axis(fig[1, 4], title = "depth-integrated TKE, eight", xlabel = L"\lambda", ylabel = L"\phi")
# hm = heatmap!(ax, λ8, ϕ8, log10.(interior(eight[:E][1], :, :, 65) .+ 1e-20), colorrange = (-5, 0), colormap = :magma, interpolate = true)

# """
# Calculate energy (Ê) and enstrophy (Ω̂) spectra on a transect in the channel
# """

# using WenoNeverworld.Diagnostics: average_spectra, hann_window

# Êcenter = average_spectra(Efluccenter, Colon(), 80:81; k = 65)
# Êweno   = average_spectra(Eflucweno,   Colon(), 80:81; k = 65)

# Ω̂center = average_spectra(Efluccenter, Colon(), 117:118; k = 65, windowing = hann_window)
# Ω̂weno   = average_spectra(Eflucweno,   Colon(), 117:118; k = 65, windowing = hann_window)

# @show ϕ[200]

# NxE = length(xnodes(weno[:v][1]))
# Nxζ = length(xnodes(weno[:u][1]))

# fx = Êcenter.freq[2:end]
# fxE = fftfreq(NxE, 1 / Δxᶜᶜᶜ(1, 1, 1, weno[:v].grid))[1:Int(NxE ÷ 2)] .* 1e3 # in 1/km
# fxζ = fftfreq(Nxζ, 1 / Δxᶠᶠᶜ(1, 1, 1, weno[:u].grid))[1:Int(Nxζ ÷ 2)] .* 1e3 # in 1/km

# fig = Figure(resolution = (800, 350))
# ax = Axis(fig[1, 1], title = "energy spectra at 50ᵒ S", yscale = log10, xscale = log10,
#                             xlabel = "Wavenumber 1/km", ylabel = L"E_k")
# lines!(ax, fxE[2:end], Êcenter.spec[2:end], color = :red,   label = "center momentum scheme")
# lines!(ax, fxE[2:end], Êweno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
# lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
# lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3),  color = :grey, linewidth = 2, linestyle = :dashdot)

# ax = Axis(fig[1, 2], title = "energy spectra at 10ᵒ S", yscale = log10, xscale = log10,
#                             xlabel = "Wavenumber 1/km", ylabel = L"E_k")
# lines!(ax, fxE[2:end], Ω̂center.spec[2:end], color = :red,   label = "center momentum scheme")
# lines!(ax, fxE[2:end], Ω̂weno.spec[2:end],   color = :blue,  label = "WENO momentum scheme")
# lines!(ax, fxE[10:end-51],   0.0006 .* fx[10:end-50].^(-5/3), color = :grey, linewidth = 2, linestyle = :dashdot)
# lines!(ax, fxE[end-100:end], 0.0055 .* fx[end-100:end].^(-3), color = :grey, linewidth = 2, linestyle = :dashdot)

# display(fig)
# CairoMakie.save("spectra2.png", fig)
# CairoMakie.save("spectra2.eps", fig)

# """ 
# Calculate MOC
# """

# ψ̄intc = WenoNeverworld.Diagnostics.calculate_residual_MOC(center_avg[:v], center_avg[:b])
# ψ̄intw = WenoNeverworld.Diagnostics.calculate_residual_MOC(weno_avg[:v],   weno_avg[:b])

# κeffc = WenoNeverworld.Diagnostics.calculate_κeff(center_avg[:b], ψ̄intc)
# κeffw = WenoNeverworld.Diagnostics.calculate_κeff(weno_avg[:b],   ψ̄intw)

# ϕ = ynodes(Face, center_avg[:b].grid)

# levels = (-0.2:6e-2:0.1)

# bwmean = time_average(weno_avg[:b]) 
# bcmean = time_average(center_avg[:b])

# Nx, Ny, Nz = size(bcmean)

# bc_surf = interior(bcmean, 15:280-15, :, 69)
# bw_surf = interior(bcmean, 15:280-15, :, 69)

# Bwmean = mean(bwmean, dims = 1)
# Bcmean = mean(bcmean, dims = 1)

# blevels = collect(0.0:0.001:0.06)

# using Statistics: quantile

# function name_contours!(ax,cplot,value)
#     beginnings = Point2f0[]
#     # First plot in contour is the line plot, first arguments are the points of the contour
#     segments = cplot.plots[1][1][]
#     #@info segments[1]
#     for (i, p) in enumerate(segments)
#         # the segments are separated by NaN, which signals that a new contour starts
#         if isnan(p)
#             push!(beginnings, segments[i-1])
#         end
#     end
#     sc = scatter!(ax, beginnings, markersize=30, color=(:white, 0.1), strokecolor=:white)
#     translate!(sc, 0, 0, 1)
#     # Reshuffle the plot order, so that the scatter plot gets drawn before the line plot
#     delete!(ax, sc)
#     delete!(ax, cplot)
#     push!(ax.scene, sc)
#     push!(ax.scene, cplot)
#     anno = text!(ax, [(string(value), p) for (i, p) in enumerate(beginnings)], 
#                        align=(:center, :center), textsize=10)

#     # move, so that text is in front
#     translate!(anno, 0, 0, 2) 
# end

# bwmin = zeros(Ny)
# bcmin = zeros(Ny)
# for j in 1:Ny
#     bwmin[j] = minimum(bwmean, dims = 1)[1, j, Nz]
#     bcmin[j] = minimum(bcmean, dims = 1)[1, j, Nz]
# end

# ϕc = ynodes(Bcmean)

# fig = Figure(resolution = (1500, 300))
# ax = Axis(fig[1, 1], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")
# hm = heatmap!(ax, ϕ, blevels, ψ̄intw ./ 1e6; colorrange = (-0.1, 0.1), colormap = :balance)
# lines!(ax, ϕc, interior(Bwmean, 1, :, 69), color = :grey, linewidth = 2)
# lines!(ax, ϕc, bwmin, color = :black, linewidth = 2)

# ax = Axis(fig[1, 2], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")
# heatmap!(ax, ϕ, blevels, ψ̄intc ./ 1e6; colorrange = (-0.1, 0.1), colormap = :balance)
# cb = Colorbar(fig[1, 3], hm)
# lines!(ax, ϕc, interior(Bcmean, 1, :, 69), color = :grey, linewidth = 2)
# lines!(ax, ϕc, bcmin, color = :black, linewidth = 2)

# ax = Axis(fig[1, 4], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")
# hm = heatmap!(ax, ϕ, blevels, (abs.(ψ̄intc) .- abs.(ψ̄intw)) ./ 1e6 ; colorrange = (-0.01, 0.01), colormap = :balance)
# cb = Colorbar(fig[1, 5], hm)
# CairoMakie.save("circulation2.png", fig)


# fig = Figure(resolution = (1500, 300))
# ax = Axis(fig[1, 1], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")
# lines!(ax, ϕc, interior(Bwmean, 1, :, 69), color = :grey, linewidth = 2)
# lines!(ax, ϕc, bwmin, color = :black, linewidth = 2)

# ax = Axis(fig[1, 2], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")
# # cb = Colorbar(fig[1, 3], hm)
# lines!(ax, ϕc, interior(Bcmean, 1, :, 69), color = :grey, linewidth = 2)
# lines!(ax, ϕc, bcmin, color = :black, linewidth = 2)

# ax = Axis(fig[1, 4], title = L"residual overturning circulation \psi Sv", ylabel = L"b ms^{-2}", xlabel = L"\phi")

# # cb = Colorbar(fig[1, 5], hm)
# CairoMakie.save("circulation2.eps", fig)

