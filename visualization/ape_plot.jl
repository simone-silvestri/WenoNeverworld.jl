using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using CairoMakie, SixelTerm
using LaTeXStrings
using JLD2

variables = ("b")
stride = 20
#=
prefix_simulation_half1 = "weno_half_original_ch"
dir_half1 = "/storage4/WenoNeverworldData/half_degree/"
fields_half1 = all_fieldtimeseries(prefix_simulation_half1, dir_half1; variables, checkpointer = true);
APE_half = Diagnostics.integral_available_potential_energy(fields_half1[:b]; stride)
filename_half = "ape_half_new.jld2"
save(filename_half, Dict("APE_half" => APE_half))

prefix_simulation_half2 = "weno_half_larger_diffusivity_ch"
dir_half2 = "/storage4/WenoNeverworldData/half_degree/"
fields_half2 = all_fieldtimeseries(prefix_simulation_half2, dir_half2; variables, checkpointer = true);
APE_half_diff = Diagnostics.integral_available_potential_energy(fields_half2[:b]; stride)
filename_half_diff = "ape_half_diff.jld2"
save(filename_half_diff, Dict("APE_half_diff" => APE_half_diff))

prefix_simulation_half3 = "weno_half_larger_viscosity_ch"
dir_half3 = "/storage4/WenoNeverworldData/half_degree/"
fields_half3 = all_fieldtimeseries(prefix_simulation_half3, dir_half3; variables, checkpointer = true);
APE_half_visc = Diagnostics.integral_available_potential_energy(fields_half3[:b]; stride)
filename_half_visc = "ape_half_visc.jld2"
save(filename_half_visc, Dict("APE_half_visc" => APE_half_visc))

prefix_simulation_half4 = "weno_half_variable_diff_ch"
dir_half4 = "/storage4/WenoNeverworldData/half_degree/"
fields_half4 = all_fieldtimeseries(prefix_simulation_half4, dir_half4; variables, checkpointer = true);
APE_half_vary = Diagnostics.integral_available_potential_energy(fields_half4[:b]; stride)
filename_vary = "ape_half_vary.jld2"
save(filename_vary, Dict("APE_half_vary" => APE_half_vary))


prefix_simulation_half5 = "weno_quarter_ch"
dir_half5 = "/storage4/WenoNeverworldData/quarter_degree/"
fields_half5 = all_fieldtimeseries(prefix_simulation_half5, dir_half5; variables, checkpointer = true);
APE_fourth2 = Diagnostics.integral_available_potential_energy(fields_half5[:b]; stride)
filename5 = "ape_fourth2.jld2"
save(filename5, Dict("APE_fourth2" => APE_fourth2)) #run with stride = 1

prefix_simulation_half6 = "weno_quarter_original_ch"
dir_half6 = "/storage4/WenoNeverworldData/quarter_degree/"
fields_half6 = all_fieldtimeseries(prefix_simulation_half6, dir_half6; variables, checkpointer = true);
APE_fourth_new = Diagnostics.integral_available_potential_energy(fields_half6[:b]; stride)
filename6 = "ape_fourth_new.jld2"
save(filename6, Dict("APE_fourth_new" => APE_fourth_new))  #run with stride = 1

=#
prefix_simulation_half7 = "weno_half_extreme_diffusivity_ch"
dir_half7 = "/storage4/WenoNeverworldData/half_degree/"
fields_half7 = all_fieldtimeseries(prefix_simulation_half7, dir_half7; variables, checkpointer = true);
APE_half_ex = Diagnostics.integral_available_potential_energy(fields_half7[:b]; stride)
filename7 = "ape_half_extreme.jld2"
save(filename7, Dict("APE_half_ex" => APE_half_ex))

prefix_simulation_half8 = "weno_half_ch"
dir_half8 = "/storage4/WenoNeverworldData/half_degree_new/"
fields_half8 = all_fieldtimeseries(prefix_simulation_half8, dir_half8; variables, checkpointer = true);
APE_half1 = Diagnostics.integral_available_potential_energy(fields_half8[:b]; stride)
filename8 = "ape_half.jld2"
save(filename8, Dict("APE_half1" => APE_half1))

fig1 = Figure(resolution = (1200, 800))
ax = Axis(fig1[1, 1], xlabel="t", xlabelsize = 20, xticklabelsize = 20, ylabel=L"[m^5/s^2]", ylabelsize = 20,title="APE 1/2 degree", yticklabelsize = 20, titlesize=25 )

#hfile1 = jldopen("ape_half_new.jld2", "r")
#hfile2 = jldopen("ape_half_diff.jld2", "r")
#hfile3 = jldopen("ape_half_visc.jld2", "r")
#hfile4 = jldopen("ape_half_vary.jld2", "r")
#hfile5 = jldopen("ape_fourth2.jld2", "r")
#hfile6 = jldopen("ape_fourth_new.jld2", "r")
hfile7 = jldopen("ape_half_extreme.jld2", "r")
hfile8 = jldopen("ape_half.jld2", "r")


#keys(hfile)
#ape_data1 = hfile1["APE_half"]          #original setup
#ape_data2 = hfile2["APE_half_diff"]     #larger diffusivity
#ape_data3 = hfile3["APE_half_visc"]     #larger viscosity and diffusivity 
#ape_data4 = hfile4["APE_half_vary"]     #variable diffusivity profile
#ape_data5 = hfile5["APE_fourth2"]        #initialized from scratch on satori
#ape_data6 = hfile6["APE_fourth_new"]     #interpolated from 1/2 degree
ape_data7 = hfile7["APE_half_ex"]
ape_data8 = hfile8["APE_half1"] 
#=
ape_data1_norm = ape_data1 / maximum(ape_data1)
ape_data2_norm = ape_data2 / maximum(ape_data2)
ape_data3_norm = ape_data3 / maximum(ape_data3)
ape_data4_norm = ape_data4 / maximum(ape_data4)
ape_data5_norm = ape_data5 / maximum(ape_data5)
ape_data6_norm = ape_data6 / maximum(ape_data6)
=#
#l1 = lines!(ax, ape_data1, color=:blue)
#l2 = lines!(ax, ape_data2, color=:red)
#l3 = lines!(ax, ape_data3, color=:green)
#l4 = lines!(ax, ape_data4, color=:purple)
#ape_data5_norm = ape_data5 / maximum(ape_data5)
#ape_data6_norm = ape_data6 / maximum(ape_data6)
#l5 = lines!(ax, ape_data5, color=:blue)
l7 = lines!(ax, ape_data7, color=:green)
l7 = lines!(ax, ape_data7, color=:red)

leg = Legend(fig1[1, 2],
    [l7, l8], ["1/2 κ=1e-3", "1/2 κ=3e-5 new"], position = :right, labelsize = 15)
    #["1/2 orginal", "1/2 larger diff", "1/2 larger visc + diff", "1/2 variable diff", "1/4 original", "1/4 larger diff"], position = :right, labelsize = 15)

display(fig1)
save("ape_compare_new.png", fig1)


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

