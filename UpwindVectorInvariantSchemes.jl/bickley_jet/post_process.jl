import Pkg
Pkg.activate("/Users/simonesilvestri/development/GlobalShenanigans.jl/")
using GlobalShenanigans
using Oceananigans
using Oceananigans.Grids: xnodes, ynodes
using Oceananigans.Fields: @compute, set!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom
using GLMakie
using JLD2

include("utils.jl")

folder = "../"

# Select cases
advection = [:WENOVel, 
             :WENOVort, 
             :WENO, 
             :WENODiv, 
             :WENODiv2, 
             :WENODiv3,
             :WENO,
             :WENO,
             :WENO,
             :WENO] 

advnames  = [:w1_, 
             :w2_, 
             :w3_, 
             :v1_, 
             :wn1_, 
             :wn2_]


varnames  = [:u, :v, :ζ, :c]
# sizes     = [:(64), :(128), :(256), :(1024)]

sizes     = [:(64)]
condition = nothing

for (adv, scheme) in zip(advnames, advection)
    for Nh in sizes
        file = Symbol(:file_, adv, Nh)
        path_str = folder * "bickley_jet_Nh_" * string(Nh) * "_" * string(scheme) * ".jld2"

        @info "Loading data file $path_str"
        @eval begin
            $file = jldopen($path_str) 
            iterations = parse.(Int, keys($file["timeseries/t"]))           
        end
        grid = RectilinearGrid(size=(Nh, Nh), x = (-2π, 2π), y=(-2π, 2π), topology = (Periodic, Periodic, Flat))
        
        if folder == "immersed"
            grid = ImmersedBoundaryGrid(grid, GridFittedBottom((x, y) -> Int(!(immersed_condition(x, y)))))
        end
        utmp = XFaceField(grid)
        vtmp = YFaceField(grid)
        ctmp = CenterField(grid)
        ζtmp = Field((Face, Face, Center), grid)

        @info "Evaluating u, v, c and ζ"
        # Defining the variables
        for (name, field) in zip(varnames, [utmp, vtmp, ctmp, ζtmp])
            var = Symbol(name, :_, adv, Nh)
            key = string(name)
            spec_r = Symbol(:spec_, name, :_, adv, Nh, :_r)
            spec_x = Symbol(:spec_, name, :_, adv, Nh, :_x)
            spec_y = Symbol(:spec_, name, :_, adv, Nh, :_y)
            @eval begin
                $var  = Vector(undef, length(iterations))
                $spec_r = Vector(undef, length(iterations))
                $spec_x = Vector(undef, length(iterations))
                $spec_y = Vector(undef, length(iterations))
                for (idx, i) in enumerate(iterations)
                    $var[idx] = deepcopy(set!($field, $file["timeseries/" * $key * "/" * string(i)][:, :, 1]))
                    $spec_r[idx] = power_spectrum_1d_isotropic(interior($var[idx])[:, :, 1], xnodes($var[idx]), ynodes($var[idx]); condition)
                    $spec_x[idx] = power_spectrum_1d_x(interior($var[idx])[:, :, 1], xnodes($var[idx]), ynodes($var[idx]); condition)
                    $spec_y[idx] = power_spectrum_1d_y(interior($var[idx])[:, :, 1], xnodes($var[idx]), ynodes($var[idx]); condition)
                end
            end
        end

        @info "Evaluating tke"
        k = Symbol(:k_, adv, Nh)
        u = Symbol(:u_, adv, Nh)
        v = Symbol(:v_, adv, Nh)
        @eval begin
            $k = Vector(undef, length(iterations))
            for (idx, i) in enumerate(iterations)
                @compute $k[idx] = Field($u[idx] ^2 + $v[idx]^2)
            end
        end
        for dir in [:_x, :_y, :_r]
            spec = Symbol(:spec_k_, adv, Nh, dir)
            sp_u = Symbol(:spec_u_, adv, Nh, dir)
            sp_v = Symbol(:spec_v_, adv, Nh, dir)
            @eval begin
                $spec = Vector(undef, length(iterations))
                for (idx, i) in enumerate(iterations)
                    $spec[idx] = Spectrum($sp_u[idx].spec .+ $sp_v[idx].spec, $sp_u[idx].freq)
                end 
            end
        end

        @info "Closing data files"
        @eval begin
            close($file)
        end
    end
end

for n in sizes
    freq  = Symbol(:f, n)
    freqr = Symbol(:f_r, n)
    orig  = Symbol(:spec_c_w1_, n, :_x)
    origr = Symbol(:spec_c_w1_, n, :_r)
    @eval begin
        $freq  = $orig[1].freq
        $freqr = $origr[1].freq
    end
end

for name in [varnames..., :k], adv in advnames, n in sizes, dir in [:_x, :_y, :_r]
    spec = Symbol(:hov_, name, :_, adv, n, dir)
    orig = Symbol(:spec_, name, :_, adv, n, dir)
    @eval begin
        $spec = zeros(Int($n/2), 101)
        for i in 1:101
            $spec[:, i] = $orig[i].spec        
        end
    end
end

####
#### Average quantities (enstr and tke)
####
using Oceananigans.AbstractOperations: volume
using Statistics: mean

for adv in advnames, n in sizes
    energy = Symbol(:energy_, adv, n)
    enstrp = Symbol(:enstrophy_, adv, n)
    k = Symbol(:k_, adv, n)
    ζ = Symbol(:ζ_, adv, n)
    
    @info "calculating mean diagnostics for $adv, $n"
    grid = RectilinearGrid(size=(n, n), x = (-2π, 2π), y=(-2π, 2π), topology = (Periodic, Periodic, Flat))
    if folder == "immersed"
        grid = ImmersedBoundaryGrid(grid, GridFittedBottom((x, y) -> Int(!(immersed_condition(x, y)))))
    end

    @eval begin
        $enstrp = Vector(undef, length($iterations))
        $energy = Vector(undef, length($iterations))
        for idx in 1:length(iterations)
            $enstrp[idx] = mean(interior($ζ[idx]) .* interior($ζ[idx]))
            $energy[idx] = mean(interior($k[idx]))
        end

        $enstrp = Float64.($enstrp)
        $energy = Float64.($energy)
    end
end


