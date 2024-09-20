using WenoNeverworld
using Oceananigans
using Oceananigans.Grids: AbstractGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

"""returns a nametuple of (u, v, w, b) from the data in file"""
function checkpoint_fields(file)
    file = jldopen(file)
    grid = file["grid"]

    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    b = CenterField(grid)
    
    Hx, Hy, Hz = halo_size(grid)
    for (var, name) in zip((u, v, w, b), ("u", "v", "w", "b"))
        set!(var, file[name * "/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz])
        fill_halo_regions!(var)
    end

    return (; u, v, w, b)
end

assumed_location(var) = var == "u" ? (Face, Center, Center) : 
                        var == "v" ? (Center, Face, Center) : 
                        var == "w" ? (Center, Center, Face) : 
                        (Center, Center, Center)

remove_last_character(s) = s[1:end-1]

"""
    all_fieldtimeseries(filename, dir = nothing; variables = ("u", "v", "w", "b"), checkpointer = false, number_files = nothing)

Load and return a dictionary of field time series data.

# Arguments
- `filename`: The name of the file containing the field data.
- `dir`: The directory where the field data files are located. Defaults to "./".
- `variables`: A tuple of variable names to load. Defaults to `("u", "v", "w", "b")`.
- `checkpointer`: A boolean indicating whether to read checkpointers or time series. Defaults to `false`.
- `number_files`: The number of files to load. Defaults to `nothing`.

# Returns
A dictionary where the keys are symbols representing the variable names and the values are `FieldTimeSeries` objects.
"""
function all_fieldtimeseries(filename, dir = "./"; 
                             variables = ("u", "v", "w", "b"), 
                             checkpointer = false,
                             number_files = nothing)

    fields = Dict()

    if !(checkpointer)
        for var in variables
            fields[Symbol(var)] = FieldTimeSeries(dir * filename, var; backend=OnDisk(), architecture=CPU())
        end
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

        times = zeros(length(numbers))

        for (idx, file) in enumerate(myfiles)
            @info "time from index $idx" file
            times[idx] = jldopen(dir * file * "2")["clock"].time
        end

        @info "loading iterations" numbers
        grid = jldopen(dir * myfiles[1] * "2")["grid"]

	# Fix the grid if we are reconstructing
	bfield = jldopen(dir * myfiles[1] * "2")["b/data"]
	Hx = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.Hx
	Hy = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.Hy
	Hz = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.Hz
	Ny = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.Ny
	Nz = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.Nz
        φF = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.φᵃᶠᵃ[1:Ny+1]
        zF = jldopen(dir * myfiles[1] * "2")["grid"].underlying_grid.zᵃᵃᶠ[1:Nz+1]

        resolution = (φF[end] - φF[1]) / Ny 
        
	grid = grid isa AbstractGrid ? grid : NeverworldGrid(resolution; z_faces = zF)

        for var in variables
            field = FieldTimeSeries{assumed_location(var)...}(grid, times)
            for (idx, file) in enumerate(myfiles)
                @info "index $idx" file
                concrete_var = jldopen(dir * file * "2")[var * "/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
		interior(field[idx]) .= concrete_var
	    end

            fields[Symbol(var)] = field
        end
    end

    return fields
end

"""limit the timeseries to `times`"""
function limit_timeseries!(fields::Dict, times)
    new_fields = Dict()

    for (key, field) in fields
        new_fields[key] = limit_timeseries!(field, times)
    end

    return new_fields
end

"""limit the timeseries to `times`"""
function limit_timeseries!(field::FieldTimeSeries, times)

    loc = location(field)
    new_field = FieldTimeSeries{loc...}(field.grid, times)

    for (idx, time) in enumerate(field.times)
        id2 = findfirst(isequal(time), times)
        if !isnothing(id2)
            set!(new_field[id2], field[idx])
        end
    end

    return new_field
end

"""saves a new file with name `new_file_name` with the last `limit_to` timeseries"""
function reduce_output_size!(old_file_name, new_file_name; limit_to = 20, variables = ("u", "v", "w", "b"))

    var_dict = all_fieldtimeseries(old_file_name; variables)
    times = var_dict[Symbol(variables[1])].times[end - limit_to:end]
    var_dict = limit_timeseries!(var_dict, times)

    jldsave(new_file_name, vars = var_dict)
end

"""adds the kinetic energy to a timeseries of values averaged in time. The latter must contains u2 and v2"""
function add_kinetic_energy_to_averaged_timeseries!(fields::Dict)

    E = FieldTimeSeries{Center, Center, Center}(fields[:u].grid, fields[:u].times[iterations])

    for i in 1:E.times
        u2 = fields[:u2][t]
        v2 = fields[:v2][t]
        set!(E[i], compute!(Field(@at (Center, Center, Center) 0.5 * (u2 + v2))))
    end

    fields[:E] = E

    return nothing
end

"""adds the kinetic energy and vertical vorticity to an instantaneous timeseries"""
function add_kinetic_energy_and_vorticity_to_timeseries!(fields::Dict)

    ζ =     FieldTimeSeries{Face, Face, Center}(fields[:u].grid, fields[:u].times)
    E = FieldTimeSeries{Center, Center, Center}(fields[:u].grid, fields[:u].times)

    for t in 1:length(E.times)
        set!(ζ[t], VerticalVorticity(fields, t))
        set!(E[t], KineticEnergy(fields, t))
    end

    fields[:E] = E
    fields[:ζ] = ζ

    return nothing
end
