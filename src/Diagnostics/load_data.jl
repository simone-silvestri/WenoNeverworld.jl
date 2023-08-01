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

"""
    all_fieldtimeseries(filename, dir = nothing; variables = ("u", "v", "w", "b"), checkpointer = false)

returns a dictionary containing a `FieldTimeSeries` for each variable in `variables`.
If `checkpointer == true` it loads the data from all the checkpoint files contained in the directory `dir`
"""
function all_fieldtimeseries(filename, dir = nothing; variables = ("u", "v", "w", "b"), checkpointer = false)

    fields = Dict()

    if !(checkpointer)
        for var in variables
            fields[Symbol(var)] = FieldTimeSeries(filename, var; backend=OnDisk(), architecture=CPU())
        end
    else
        files = readdir(dir)
        files = filter((x) -> length(x) >= length(filename), files)
        myfiles = filter((x) -> x[1:length(filename)] == filename, files)
        numbers = parse.(Int, filter.(isdigit, myfiles))
        times   = numbers
        
        @info "loading times" times
        grid = jldopen(dir * myfiles[1])["grid"] 
        perm = sortperm(numbers)
        myfiles = myfiles[perm]
        for var in variables
            field = FieldTimeSeries{assumed_location(var)...}(grid, times)
            for (idx, file) in enumerate(myfiles)
                @info "index $idx" file
                concrete_var = jldopen(dir * file)[var * "/data"]
                set!(field[idx], concrete_var)
            end

            fields[Symbol(var)] = field
        end
    end

    return fields
end

"""like the function above but extract the surface"""
function all_fieldtimeseries_surface(filename, dir = nothing; variables = ("u", "v", "w", "b"))

    fields = Dict()
    files = readdir(dir)
    files = filter((x) -> length(x) >= length(filename), files)
    myfiles = filter((x) -> x[1:length(filename)] == filename, files)
    numbers = parse.(Int, filter.(isdigit, myfiles))
    
    @info "loading times" times
    grid = jldopen(dir * myfiles[1])["grid"] 
    Hz   = halo_size(grid)[3]
    Nz   = size(grid, 3)
    perm = sortperm(numbers)
    myfiles = myfiles[perm]
    for var in variables
        field = Field[]
        for (idx, file) in enumerate(myfiles)
            @info "index $idx" file
            try
                concrete_var = jldopen(dir * file)[var * "/data"][:, :, Hz + Nz]
                push!(field, concrete_var)
            catch err
                @warn err
        end
        fields[Symbol(var)] = field
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
        set!(ζ[t], VerticalVorticityField(fields, t))
        set!(E[t], KineticEnergyField(fields, t))
    end

    fields[:E] = E
    fields[:ζ] = ζ

    return nothing
end