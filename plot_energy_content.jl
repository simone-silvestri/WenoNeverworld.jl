using WenoNeverworld
using WenoNeverworld.Diagnostics
using GLMakie
using LaTeXStrings
using JLD2

variables = ("b")
stride = 20
#=
prefix_simulation_eighth = "weno_eighth_ch"
dir_eighth = "/storage3/WenoNeverworldData/"
fields_eighth = all_fieldtimeseries(prefix_simulation_eighth, dir_eighth; variables, checkpointer = true);
APE_eighth = Diagnostics.integral_available_potential_energy(fields_eighth[:b]; stride)

prefix_simulation_fourth = "weno_fourth_ch"
dir_fourth = "/storage2/WenoNeverworldData/"
fields_fourth = all_fieldtimeseries(prefix_simulation_fourth, dir_fourth; variables, checkpointer = true);
APE_fourth = Diagnostics.integral_available_potential_energy(fields_fourth[:b]; stride)

filename_fourth = "ape_fourth.jld2"
save(filename_fourth, APE_fourth)
=#
#stride = 8
prefix_simulation_half = "weno_half_ch"
dir_half = "/storage2/WenoNeverworldData/"
fields_half = all_fieldtimeseries(prefix_simulation_half, dir_half; variables, checkpointer = true);
APE_half = Diagnostics.integral_available_potential_energy(fields_half[:b]; stride)

filename_half = "ape_half.jld2"
save(filename_half, Dict("APE_half" => APE_half))

fig = Figure()
ax = Axis(fig[1, 1], xlabel="t", xlabelsize = 20, xticklabelsize = 20, ylabel=L"[m^5/s^2]", ylabelsize = 20,title="Available Potential Energy",  yticklabelsize = 20, titlesize=25, yticks=range(-6.7*10^23,-6.69*10^23, length =2 ))

hfile = jldopen("ape_half.jld2", "r")
#keys(hfile)
ape_data = hfile["APE_half"]
lines!(ax, ape_data, color=:blue)
#lines!(ax, APE_fourth, color=:red)
#lines!(ax, APE_eighth, color=:green)

display(fig)
save("plotting/ape_test.png", fig)


#integrated_heat_content = Diagnostics.heat_content(fields[:b]; stride)

#kinetic_energy = Diagnostics.integral_kinetic_energy(fields[:u], fields[:v])
# transport = Diagnostics.ACC_transport(fields[:u])
#i = length(fields[:u]) - 10
#KE = Diagnostics.KineticEnergy(fields, i)
#vert_vorticity = Diagnostics.VerticalVorticity(fields, i)
#heatmap(interior)?
#rossby_radius also has one

#lines((APE))
#GLMakie.lines(APE)
#GLMakie.lines(integrated_heat_content)
#GLMakie.lines(kinetic_energys)
