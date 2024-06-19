using WenoNeverworld, Oceananigans
using WenoNeverworld.Diagnostics
using WenoNeverworld.Diagnostics: Spectrum, average_spectra, VolumeField
using WenoNeverworld.Diagnostics: propagate, time_average, KineticEnergyField, PotentialEnergyField
using CairoMakie
using FFTW
using Statistics

# Define function to compute kinetic energy spectrum

function KE_integral_function(fields) 
    # Compute kinetic energy field
   
    u = fields[:u]
    v = fields[:v]

    Vol = VolumeField(u.grid)

    integrated_KE(u, v, V) = sum(interior(KineticEnergyField((; u, v))) .* interior(V), dims = (1, 3))
    
    KE_series = propagate(u, v, Vol; func = integrated_KE)
 
    @info "time averaging"
    
    KE = mean(KE_series)

    return KE
end

function APE_integral_function(fields)

    b   = fields[:b]
    Vol = VolumeField(b.grid)

    integrated_APE(b, V) = sum(interior(PotentialEnergyField((; b))) .* interior(V); dims = (1, 3))

    APE_series = propagate(b, Vol; func = integrated_APE)

    @info "time averaging"

    APE = mean(APE_series)

    return APE 
end

# Load data for different grid resolutions
prefixes = ["weno_half_ch", "weno_fourth_ch", "weno_eighth_ch", "weno_sixteenth_ch", "weno_thirtytwo_comp"]
dirs = ["/storage2/WenoNeverworldData/", "/storage2/WenoNeverworldData/", "/storage3/WenoNeverworldData/", "/storage3/WenoNeverworldData/", "/storage4/WenoNeverworldData/"]
resolutions = ["1/2", "1/4", "1/8", "1/16", "1/32"]
colors = [:red3, :darkorange, :green, :purple, :navy]

KE_integral  = []
APE_integral = []

####
#### Half degree
####

i  = 1 

fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
KE_tmp = KE_integral_function(fields)
push!(KE_integral, KE_tmp)

GC.gc(true)

fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("b", ), checkpointer = true, number_files = 1)
APE_tmp = APE_integral_function(fields)
push!(APE_integral, APE_tmp)

GC.gc(true)

####
#### Quarter degree
####

i  = 2

fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
KE_tmp = KE_integral_function(fields)
push!(KE_integral, KE_tmp)

GC.gc(true)

fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("b", ), checkpointer = true, number_files = 1)
APE_tmp = APE_integral_function(fields)
push!(APE_integral, APE_tmp)

GC.gc(true)

# ####
# #### Eighth degree
# ####
# 
i  = 3 
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
KE_tmp = KE_integral_function(fields)
push!(KE_integral, KE_tmp)
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("b", ), checkpointer = true, number_files = 1)
APE_tmp = APE_integral_function(fields)
push!(APE_integral, APE_tmp)
# 
GC.gc(true)
# 
# ####
# #### Sixteenth degree
# ####
# 
i  = 4 
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
KE_tmp = KE_integral_function(fields)
push!(KE_integral, KE_tmp)
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("b", ), checkpointer = true, number_files = 1)
APE_tmp = APE_integral_function(fields)
push!(APE_integral, APE_tmp)
# 
GC.gc(true)
# 
i  = 5
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("u", "v"), checkpointer = true, number_files = 1)
KE_tmp = KE_integral_function(fields)
push!(KE_integral, KE_tmp)
# 
fields = all_fieldtimeseries(prefixes[i], dirs[i]; variables = ("b", ), checkpointer = true, number_files = 1)
APE_tmp = APE_integral_function(fields)
push!(APE_integral, APE_tmp)
# 
GC.gc(true)
# 
# 
# #########
# ##### Figure creation
# ########


fig = Figure(resolution = (1000, 800))
ax = Axis(fig[1, 1], yticklabelsize = 25, xticklabelsize = 25)

lines!(ax, KE_integral[1][1, :, 1], color = colors[1])
lines!(ax, KE_integral[2][1, :, 1], color = colors[2])
lines!(ax, KE_integral[3][1, :, 1], color = colors[3])
lines!(ax, KE_integral[4][1, :, 1], color = colors[4])
lines!(ax, KE_integral[5][1, :, 1], color = colors[5])


# display(fig)
save("ke_int2.png", fig)


fig2 = Figure(resolution = (1000, 800))
ax2 = Axis(fig2[1, 1], yticklabelsize = 25, xticklabelsize = 25)
lines!(ax2, APE_integral[1][1, :, 1], color = colors[1])
lines!(ax2, APE_integral[2][1, :, 1], color = colors[2])
lines!(ax2, APE_integral[3][1, :, 1], color = colors[3])
lines!(ax2, APE_integral[4][1, :, 1], color = colors[4])
lines!(ax2, APE_integral[5][1, :, 1], color = colors[5])

display(fig2)
save("ape_int2.png", fig2)
