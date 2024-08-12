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

prefix_simulation_half2 = "weno_half_original_ch"
dir_half2 = "/storage4/WenoNeverworldData/half_degree/"
fields_half2 = all_fieldtimeseries(prefix_simulation_half2, dir_half2; variables, checkpointer = true);
APE_half = Diagnostics.integral_available_potential_energy(fields_half2[:b]; stride)
filename_half = "ape_half_new.jld2"
save(filename_half, Dict("APE_half" => APE_half))
=#
prefix_simulation_half1 = "weno_half_larger_diffusivity_ch"
dir_half1 = "/storage4/WenoNeverworldData/half_degree/"
fields_half1 = all_fieldtimeseries(prefix_simulation_half1, dir_half1; variables, checkpointer = true);
APE_half_diff = Diagnostics.integral_available_potential_energy(fields_half1[:b]; stride)
filename_half_warm = "ape_half_diff.jld2"
save(filename_half_warm, Dict("APE_half_diff" => APE_half_diff))

prefix_simulation_half2 = "weno_half_larger_viscosity_ch"
dir_half2 = "/storage4/WenoNeverworldData/half_degree/"
fields_half2 = all_fieldtimeseries(prefix_simulation_half2, dir_half2; variables, checkpointer = true);
APE_half_visc = Diagnostics.integral_available_potential_energy(fields_half2[:b]; stride)
filename_half = "ape_half_visc.jld2"
save(filename_half, Dict("APE_half_visc" => APE_half_visc))

prefix_simulation_half4 = "weno_half_variable_diff_ch"
dir_half4 = "/storage4/WenoNeverworldData/half_degree/"
fields_half4 = all_fieldtimeseries(prefix_simulation_half4, dir_half4; variables, checkpointer = true);
APE_half_vary = Diagnostics.integral_available_potential_energy(fields_half4[:b]; stride)
filename_vary = "ape_half_vary.jld2"
save(filename_vary, Dict("APE_half_vary" => APE_half_vary))


fig2 = Figure(resolution = (1200, 800))
ax = Axis(fig2[1, 1], xlabel="t", xlabelsize = 20, xticklabelsize = 20, ylabel=L"[m^5/s^2]", ylabelsize = 20,title="Available Potential Energy",  yticklabelsize = 20, titlesize=25, yticks=range(-6.7*10^23,-6.69*10^23, length =2 ))

hfile1 = jldopen("ape_half_new.jld2", "r")
hfile2 = jldopen("ape_half_diff.jld2", "r")
hfile3 = jldopen("ape_half_visc.jld2", "r")
hfile4 = jldopen("ape_half_vary.jld2", "r")

#keys(hfile)
ape_data1 = hfile1["APE_half"]          #original setup
ape_data2 = hfile2["APE_half_diff"]     #larger diffusivity
ape_data3 = hfile3["APE_half_visc"]     #larger viscosity and diffusivity 
ape_data4 = hfile4["APE_half_vary"]     #variable diffusivity profile
l1 = lines!(ax, ape_data1, color=:blue)
l2 = lines!(ax, ape_data2, color=:red)
l3 = lines!(ax, ape_data3, color=:green)
l4 = lines!(ax, ape_data4, color=:purple)

leg = Legend(fig2[1, 2],
    [l1, l2, l3, l4],
    ["orginal", "larger diff", "larger visc + diff", "variable diff"], position = :right, labelsize = 15)

display(fig2)
save("plotting/ape_all1.png", fig2)


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
