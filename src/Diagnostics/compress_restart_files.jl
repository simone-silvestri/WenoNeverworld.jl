
function check_ranges(folder, ranks; H = 7, iteration = 0)
    Ny = Vector(undef, ranks)
    jranges = Vector(undef, ranks)
    for rank in 0:ranks - 1
        var = jldopen(folder * output_prefix * "0_checkpoint_iteration$(iteration).jld2")["u/data"][H+1:end-H, H+1:end-H, H+1:end-H] 
        Ny[rank+1] = size(var, 2)
    end
    jranges[1] = UnitRange(1, Ny[1])
    for rank in 2:ranks
        jranges[rank] = UnitRange(jranges[rank-1][end]+1,jranges[rank-1][end]+Ny[rank])
    end

    return jranges
end

function compress_restart_file(resolution, ranks, iteration, folder = "../"; 
                               output_prefix = "weno_thirtytwo",
                               full_grid = NeverworldGrid(resolution))

    @info "initializing active map"
    bathymetry = interior(full_grid.immersed_boundary.bottom_heigth, :, :, 1)

    fields_data = Dict()
    fields_data[:underlying_grid] = full_grid
    fields_data[:bathymetry] = bathymetry
    fields_data[:clock] = jldopen(folder * output_prefix * "0_checkpoint_iteration$(iteration).jld2")["clock"]

    jranges = check_ranges(folder, ranks; output_prefix, H, iteration)

    @info "starting the compression"
    for var in (:u, :w, :v, :T, :S)
        GC.gc()

        @info "compressing variable $var"
        sizefield = var == :v ? (Nx, Ny+1, Nz) :
                    var == :w ? (Nx, Ny, Nz+1) : (Nx, Ny, Nz)

        compressed_data = zeros(Float32, sizefield...)

        for rank in 0:ranks-1
            @info "reading rank $rank"
            jrange = iranges[rank+1]
            compressed_data[:, jrange, :] .= jldopen(folder * output_prefix * "$(rank)_checkpoint_iteration$(iteration).jld2")[string(var) * "/data"][H+1:end-H, H+1:end-H, H+1:end-H]
        end

        fields_data[var] = compressed_data
    end

    compressed_η = zeros(Float32, Nx, Ny, 1)
    for rank in 0:ranks-1
        @info "reading rank $rank"

        jrange = jranges[rank+1]
        data = jldopen(folder * output_prefix * "$(rank)_checkpoint_iteration$(iteration).jld2")["η/data"]
        Hx = calc_free_surface_halo(irange, data)
        data = data[Hx+1:end-Hx, H+1:end-H, :]
        compressed_η[:, jrange, :] .= Float32.(data)
    end

    fields_data[:η] = compressed_η

    jldopen(folder * "compressed_iteration_$(iteration).jld2","w") do f
        for (key, value) in fields_data
            f[string(key)] = value
        end
    end
end

function calc_free_surface_halo(jrange, data)
    Ny = size(data, 2)
    ny = length(jrange)
    return Int((Ny - ny) ÷ 2)
end

const regex = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";

function compress_all_restarts(resolution, ranks, dir; remove_restart = false)
    files = readdir(dir)
    files = filter(x -> length(x) > 30, files)
    files = filter(x -> x[1:26] == "RealisticOcean_checkpoint_", files)
    iterations = Int[]
    for file in files
        file   = file[1:end-5]
        string = ""
        i = length(file)
        while occursin(regex, "$(file[i])")
            string = file[i] * string
            i -= 1
        end
        push!(iterations, parse(Int, string))
    end

    iterations = unique(iterations)
    iterations = sort(iterations)
    for iter in iterations
        @info "compressing iteration $iter"
        compress_restart_file(resolution, ranks, iter, dir)

        if remove_restart
            @info "removing iteration $iter"

            for rank in 0:ranks-1
                to_remove = "RealisticOcean_checkpoint_$(rank)_iteration$(iter).jld2"
                cmd = `rm $to_remove`
                run(cmd)
            end
        end
    end

    return nothing
end