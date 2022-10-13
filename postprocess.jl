using GLMakie
using JLD2

include("initial_conditions.jl")
include("neverworld_utils.jl")

using Oceananigans.Fields: @compute

# function assumed_location(a::String) 
#     if a[1] == 'u' || a[1] == 'k'
#         return :((Face, Center, Center))
#     elseif  a[1] == 'v' 
#         return :((Center, Face, Center)) 
#     elseif  a[1] == 'w' 
#         return :((Center, Center, Face)) 
#     elseif  a[1] == 'b' 
#         return :((Center, Center, Center)) 
#     elseif  a[1] == 'η' 
#         return :((Center, Center, Nothing)) 
#     elseif a[1] == 'ζ'
#         return :((Face, Face, Center)) 
#     elseif a[1] == 't'
#         return :((Nothing, Nothing, Nothing)) 
#     else
#         return :((Center, Center, Center)) 
#     end
# end

# filename1 = "final_quarter/neverworld_quarter_averaged.jld2"
# filename2 = "final_eighth/neverworld_eighth_averaged.jld2"
# filename3 = "sixteenth-first-5-years/neverworld_sixteenth_averaged.jld2"

# for (idx, filename, res) in zip((1, 2, 3), (filename1, filename2, filename3), (1/4, 1/8, 1/16)) 

#     grid = neverworld_grid(CPU(), res; H = 3)

#     file = jldopen(filename)
#     iteration = parse(Int, keys(file["timeseries/u"])[end])

#     for key in keys(file["timeseries"])
#         dat = Symbol(key, :_data)
#         var = Symbol(key, :_, idx)
#         avg = Symbol(key, :avg, :_, idx)
#         int = Symbol(key, :int, :_, idx)
#         int1 = Symbol(key, :int1, :_, idx)

#         loc = assumed_location(key)
#         @info "reading " * key
#         @eval begin
#             $dat = $file["timeseries/$($key)/$($iteration)"]
#             if $dat isa String
#                 $dat = parse(Float64, $dat)
#             end
#             $var  = Field($loc, $grid)

#             @show size($var), size($dat)
#             Oceananigans.Fields.set!($var, $dat)
#             @compute $avg  = Field(Integral($var, dims = 1))
#             @compute $int  = Field(Integral($var, dims = 3))
#             @compute $int1 = Field(Integral($var, dims = (1, 3)))
#         end
#     end
#     close(file)
# end

# # calculating TKE
# @compute tke_1 = Field(u2_1^2 + v2_1^2 + w2_1^2)
# @compute tke_2 = Field(u2_2^2 + v2_2^2 + w2_2^2)
# @compute tke_3 = Field(u2_3^2 + v2_3^2 + w2_3^2)

# @compute tkeint_1 = Field(Integral(tke_1, dims = (1, 3)));
# @compute tkeint_2 = Field(Integral(tke_2, dims = (1, 3)));
# @compute tkeint_3 = Field(Integral(tke_3, dims = (1, 3)));

# compute the timeseries to check that we are converged

function compute_convergence(folder, sim, grid; H = 4)
    l = length("neverworld_" * sim * "_checkpoint_iteration")
    k  = []
    it = []
    for file in readdir(folder)
        if length(file) > l
            if file[1:l] == "neverworld_" * sim * "_checkpoint_iteration"
                iteration = parse(Int, file[l+1:end-5])

                f = jldopen(folder * file)

                @info "doing iteration $iteration"
                u, v = out_var(f, grid, H)

                kfield = compute!(Field(u^2 + v^2))
                push!(k, compute!(Field(Integral(kfield)))[1, 1, 1])
                push!(it, iteration)

                close(f)
            end
        end
    end
    p = sortperm([it...])
    return [k...][p], [it...][p]
end

# file_quarter   = jldopen("final_quarter/neverworld_quarter_checkpoint_iteration8213138.jld2")
# file_sixteenth = jldopen("sixteenth-first-5-years/neverworld_sixteenth_checkpoint_iteration989578.jld2")

# grid16 = neverworld_grid(CPU(), 1/16; H=4);
# grid4  = neverworld_grid(CPU(), 1/4; H=3);

# function out_var(file, grid, H)
#     ud = file["u/data"][H+1:end-H, H+1:end-H, H+1:end-H]
#     vd = file["v/data"][H+1:end-H, H+1:end-H, H+1:end-H]
    
#     u, v, w = VelocityFields(grid)
#     @show size(u), size(ud)
#     set!(u, ud)
#     set!(v, vd)

#     return (u, v)
# end

# using Oceananigans.Operators

# k4  = compute!(Field(u4^2 + v4^2))
# k16 = compute!(Field(u16^2 + v16^2))

# kint4  = Field(Integral(k4, dims = (1, 3)));
# kint16 = Field(Integral(k16, dims = (1, 3)));