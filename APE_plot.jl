using WenoNeverworld, Oceananigans
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute
using GLMakie
using LaTeXStrings
using JLD2
using HDF5

variables = ("b",)

assumed_location(var) = (Center, Center, Center)
remove_last_character(s) = s[1:end-1]

function b_fieldtimeseries(filename, dir = nothing;
                             checkpointer = false,
                             number_files = nothing)

    fields = Dict()

    if !(checkpointer)
        fields[b] = FieldTimeSeries(filename, "b"; backend=OnDisk(), architecture=CPU())
    else
        files = readdir(dir)
        files = filter((x) -> length(x) >= length(filename), files)
        myfiles = filter((x) -> x[1:length(filename)] == filename, files)
        myfiles = remove_last_character.(myfiles)
        numbers = parse.(Int, filter.(isdigit, myfiles))
        perm    = sortperm(numbers)
        numbers = numbers[perm]
        myfiles = myfiles[perm]

        if !isnothing(number_files)
            numbers = numbers[end-number_files:end]
            myfiles = myfiles[end-number_files:end]
        end

        @info "loading iterations" numbers

        # out of memory error here for 1/16
        temp = jldopen(dir * myfiles[1] * "2")

        grid = try
            jldopen(dir * myfiles[1] * "2")["grid"] 
        catch
            NeverworldGrid(jldopen(dir * myfiles[1] * "2")["resolution"])
        end
            
        field = FieldTimeSeries{assumed_location("b")...}(grid, numbers)
        for (idx, file) in enumerate(myfiles)
            @info "index $idx" file
            concrete_var = jldopen(dir * file * "2")["b" * "/data"]
            #field.times[idx] = try
                #jldopen(dir * file * "2")["clock"].time
            #catch 
                #jldopen(dir * file * "2")["clock"].time
            #end
            #set!(field[idx], concrete_var)
            field.times[idx] = jldopen(dir * file * "2")["clock"].time
            set!(field[idx], concrete_var)

        end

        fields[Symbol("b")] = field
    end

    return fields
end

# Define a function to load data and compute APE
function compute_and_save_APE(prefix, directory, variables, stride, filename)
    fields = b_fieldtimeseries(prefix, directory; checkpointer = true)
    APE = Diagnostics.integral_available_potential_energy(fields[:b]; stride)
    save(filename, Dict("APE" => APE))
end

# Define filenames for saving APE data
#filename_half = "ape_half_test.jld2"
#filename_fourth = "ape_fourth_test.jld2"
#filename_eighth = "ape_eighth_test2.jld2"
filename_sixteen = "ape_sixteen_test.jld2"


#"test" = the correct files, aka the ones run with stride = 1 and all files (for 1/2, 1/4, 1/16)
#for the 1/8, ape_eighth_test = the storage 4 files only 7 total
#ape_eighth_test2 = all 260 files, stride = 15
#ape_spectrum_test = with just 7 eighth degree files
#ape_spectrum_test2 plot = with more eighth degree files

#stride = 5
#=
# Load and save APE data for half resolution
prefix_simulation_half = "weno_half_ch"
dir_half = "/storage2/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_half, dir_half, variables, stride, filename_half)

# Load and save APE data for fourth resolution
prefix_simulation_fourth = "weno_fourth_ch"
dir_fourth = "/storage4/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_fourth, dir_fourth, variables, stride, filename_fourth)

# Load and save APE data for fourth resolution
prefix_simulation_eighth = "weno_eighth_ch"
dir_eighth = "/storage3/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_eighth, dir_eighth, variables, stride, filename_eighth)
=#
stride = 1
prefix_simulation_sixteen = "weno_sixteenth_ch"
dir_sixteen = "/storage3/WenoNeverworldData/"
compute_and_save_APE(prefix_simulation_sixteen, dir_sixteen, variables, stride, filename_sixteen)

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