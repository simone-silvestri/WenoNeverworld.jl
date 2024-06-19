using WenoNeverworld, Oceananigans
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute

using GLMakie
using LaTeXStrings
using JLD2
using HDF5

#prefix_simulation = "weno_half_ch"
#dir= "/storage2/WenoNeverworldData/"
#variables = ("u", "v","b",)
#stride = 15
#fields = all_fieldtimeseries(prefix_simulation, dir; variables, checkpointer = true);
@info "Loading data..."
path = pwd()
hfile = jldopen("/storage2/WenoNeverworldData/weno_half_checkpoint_iteration985500.jld2", "r")
keys(hfile)
oceangrid = hfile["grid"]
halo = 7

z = oceangrid.underlying_grid.zᵃᵃᶜ[1:end-halo]
Δz = oceangrid.underlying_grid.Δzᵃᵃᶜ[1:end-halo]
lat = oceangrid.underlying_grid.φᵃᶜᵃ[1:end-halo]
lon = oceangrid.underlying_grid.λᶜᵃᵃ[1:end-halo]
η = hfile["η"]["data"][halo:end-halo, halo:end-halo]
b = hfile["b"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
u = hfile["u"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = hfile["v"]["data"][halo+1:end-halo, halo+1:end-halo, halo+1:end-halo]
v = 0.5 * (v[:, 1:end-1, :] + v[:, 2:end, :])


fields = Dict(:u => u, :v => v)
i = 1  # Choose the appropriate index
#vort = VerticalVorticity(fields, i)

#i = length(fields[:u]) - 10
vort = VerticalVorticity(fields[:u], fields[:v])  #(u, v, Δz)
#vort = Diagnostics.VerticalVorticity(u, v, Δz)
##
fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/4∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:10:70, yticksize = 15, xticksize = 15)
hm = heatmap!(ax, lon, lat, vort[:,:,1], colorrange = (-4, 1), colormap = (Reverse:plasma))
#cbar1 = Colorbar(fig[1,2], hm, width = 30, ticksize = 10, ticklabelsize = 20)
display(fig)
#save("plotting/tke_slice_eighth_colorbar.png", fig)

using CairoMakie
CairoMakie.activate!()
CairoMakie.save("vort_fourth.png", fig, px_per_unit = 5)


#=
# Define a function to load data and compute APE
function compute_and_save_vort(prefix, directory, variables, stride, filename)
    fields = b_fieldtimeseries(prefix, directory; checkpointer = true)
    vort = Diagnostics.VerticalVorticity(fields[:u], fields[:v]; i)
    A#PE = Diagnostics.integral_available_potential_energy(fields[:]; stride)
    save(filename, Dict("vort" => vort))
end

# Define filenames for saving APE data
#filename_half = "vort_half_test.jld2"
filename_fourth = "vort_fourth_test.jld2"
#filename_eighth = "ape_eighth_test2.jld2"
#filename_sixteen = "ape_sixteen_test.jld2"


#"test" = the correct files, aka the ones run with stride = 1 and all files (for 1/2, 1/4, 1/16)
#for the 1/8, ape_eighth_test = the storage 4 files only 7 total
#ape_eighth_test2 = all 260 files, stride = 15
#ape_spectrum_test = with just 7 eighth degree files
#ape_spectrum_test2 plot = with more eighth degree files

#stride = 5

# Load and save APE data for half resolution
#prefix_simulation_half = "weno_half_ch"
#dir_half = "/storage2/WenoNeverworldData/"
#compute_and_save_APE(prefix_simulation_half, dir_half, variables, stride, filename_half)

# Load and save APE data for fourth resolution
prefix_simulation_fourth = "weno_fourth_ch"
dir_fourth = "/storage4/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_fourth, dir_fourth, variables, stride, filename_fourth)

# Load and save APE data for fourth resolution
prefix_simulation_eighth = "weno_eighth_ch"
dir_eighth = "/storage3/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_eighth, dir_eighth, variables, stride, filename_eighth)

stride = 1
prefix_simulation_sixteen = "weno_sixteenth_ch"
dir_sixteen = "/storage3/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_sixteen, dir_sixteen, variables, stride, filename_sixteen)
=#
#=
using JLD2
using HDF5
using GLMakie
# axis def
minval = -6.7*10^23
maxval = -6.69*10^23
# yticks=range(minval,maxval, length =10)
fig = Figure()
ax = Axis(fig[1, 1], xlabel="t", xlabelsize = 20, xticklabelsize = 20, ylabel=L"[m^5/s^2]", ylabelsize = 20,title="Available Potential Energy",  yticklabelsize = 20, titlesize=25)

hfile_1 = jldopen("ape_half_test.jld2", "r")
ape_data_half = hfile_1["APE"]

hfile_2 = jldopen("ape_fourth_test.jld2", "r")
ape_data_fourth = hfile_2["APE"]

hfile_3 = jldopen("ape_eighth_test2.jld2", "r")
ape_data_eighth = hfile_3["APE"]

hfile_4 = jldopen("ape_sixteen.jld2", "r")
ape_data_sixteen = hfile_4["APE"]

lines!(ax, ape_data_half, color=:blue)
lines!(ax, ape_data_fourth, color=:red)
lines!(ax, ape_data_eighth, color=:green)
lines!(ax, ape_data_sixteen, color=:black)

display(fig)
save("plotting/ape_spectrum_latest.png", fig)
=#

