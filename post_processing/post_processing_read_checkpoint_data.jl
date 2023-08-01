using Pkg
using Oceananigans
using Oceananigans.Grids
using JLD2
using CairoMakie
using Printf
using WenoNeverworld


include("post_processing_make_plots.jl")


function ExtractInterior1DArray(ArrayWithHalos, interior_size, halo_size, ArrayLocation, 
                                ArrayWithHalosUsingNegativeIndex = true)
    
    interior_index_lower = halo_size .+ 1
    interior_index_upper = interior_index_lower .+ interior_size .- 1
    
    if ArrayWithHalosUsingNegativeIndex
    
        if ArrayLocation == "caa"
            return ArrayWithHalos[1:interior_size[1]]
            
        elseif ArrayLocation == "faa"
            return ArrayWithHalos[1:interior_size[1]+1]
            
        elseif ArrayLocation == "aca"
            return ArrayWithHalos[1:interior_size[2]]
            
        elseif ArrayLocation == "afa"
            return ArrayWithHalos[1:interior_size[2]+1]
            
        elseif ArrayLocation == "aac"
            return ArrayWithHalos[1:interior_size[3]]
            
        elseif ArrayLocation == "aaf"
            return ArrayWithHalos[1:interior_size[3]+1]
        
        end
    
    else
    
        if ArrayLocation == "caa"
            return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1]]
            
        elseif ArrayLocation == "faa"
            return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1]+1]
            
        elseif ArrayLocation == "aca"
            return ArrayWithHalos[interior_index_lower[2]:interior_index_upper[2]]
            
        elseif ArrayLocation == "afa"
            return ArrayWithHalos[interior_index_lower[2]:interior_index_upper[2]+1]
            
        elseif ArrayLocation == "aac"
            return ArrayWithHalos[interior_index_lower[3]:interior_index_upper[3]]
            
        elseif ArrayLocation == "aaf"
            return ArrayWithHalos[interior_index_lower[3]:interior_index_upper[3]+1]
        
        end
    
    end
    
end


function ExtractInterior3DArray(ArrayWithHalos, interior_size, halo_size, ArrayLocation)
    
    interior_index_lower = halo_size .+ 1
    interior_index_upper = interior_index_lower .+ interior_size .- 1
    
    if ArrayLocation == "ccc"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1], 
                              interior_index_lower[2]:interior_index_upper[2], 
                              interior_index_lower[3]:interior_index_upper[3]]
        
    elseif ArrayLocation == "fcc"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1] + 1, 
                              interior_index_lower[2]:interior_index_upper[2], 
                              interior_index_lower[3]:interior_index_upper[3]]
    
    elseif ArrayLocation == "cfc"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1], 
                              interior_index_lower[2]:interior_index_upper[2] + 1, 
                              interior_index_lower[3]:interior_index_upper[3]]
    
    elseif ArrayLocation == "ccf"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1], 
                              interior_index_lower[2]:interior_index_upper[2], 
                              interior_index_lower[3]:interior_index_upper[3] + 1]
        
    elseif ArrayLocation == "ffc"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1] + 1, 
                              interior_index_lower[2]:interior_index_upper[2] + 1, 
                              interior_index_lower[3]:interior_index_upper[3]]
    
    elseif ArrayLocation == "fcf"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1] + 1, 
                              interior_index_lower[2]:interior_index_upper[2], 
                              interior_index_lower[3]:interior_index_upper[3] + 1]
        
    elseif ArrayLocation == "cff"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1], 
                              interior_index_lower[2]:interior_index_upper[2] + 1, 
                              interior_index_lower[3]:interior_index_upper[3] + 1]
        
    elseif ArrayLocation == "fff"
        return ArrayWithHalos[interior_index_lower[1]:interior_index_upper[1] + 1, 
                              interior_index_lower[2]:interior_index_upper[2] + 1, 
                              interior_index_lower[3]:interior_index_upper[3] + 1]
    
    end

end


function CreateFieldsFromCheckPointOutput(grid, interior_size, halo_size, filename)
    
    file = jldopen(filename)
    T_Array = file["T"]["data"]
    T_Field = Field{Center, Center, Center}(grid)
    v_Array = file["v"]["data"]
    v_Field = Field{Center, Face, Center}(grid)
    
    interior_index_lower = halo_size .+ 1
    interior_index_upper = interior_index_lower .+ interior_size .- 1

    set!(T_Field, T_Array[interior_index_lower[1]:interior_index_upper[1], 
                          interior_index_lower[2]:interior_index_upper[2], 
                          interior_index_lower[3]:interior_index_upper[3]])
    
    set!(v_Field, v_Array[interior_index_lower[1]:interior_index_upper[1], 
                          interior_index_lower[2]:interior_index_upper[2]+1, 
                          interior_index_lower[3]:interior_index_upper[3]])
    
    return T_Field, v_Field
    
end


function ComputeStreamFunctionAndPlotMeridionalOverturningCirculation_1(path, first_index, last_index, interval_index)

    arch = CPU()
    new_degree = 1
    grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))
    
    interior_size = (grid.Nx, grid.Ny, grid.Nz) 
    halo_size = (grid.Hx, grid.Hy, grid.Hz)
    
    φᵃᶠᵃ_Array_Interior = ExtractInterior1DArray(grid.φᵃᶠᵃ, interior_size, halo_size, "afa")
    φᵃᶠᵃ_Array_Interior_Plot = zeros(grid.Ny+1)
    φᵃᶠᵃ_Array_Interior_Plot[:] = φᵃᶠᵃ_Array_Interior[:]
    
    zᵃᵃᶠ_Array_Interior = ExtractInterior1DArray(grid.zᵃᵃᶠ, interior_size, halo_size, "aaf")
    zᵃᵃᶠ_Array_Interior_Plot = zeros(grid.Nz+1)
    zᵃᵃᶠ_Array_Interior_Plot[:] = zᵃᵃᶠ_Array_Interior[:]/1000
    
    ψ_Plot = zeros(grid.Ny+1, grid.Nz+1)
    n_indices = Int((last_index - first_index)/interval_index) + 1
    ψ_Mean_Plot = zeros(grid.Ny+1, grid.Nz+1)
    
    make_animation = false
    if make_animation
        ψ_Plot_TimeSeries = zeros(n_indices, grid.Ny+1, grid.Nz+1)
        title_TimeSeries = fill("", n_indices)
    end
    
    int_T_xyz_TimeSeries = zeros(n_indices)
    iTimeSeries = 0
    
    resolution = (900, 750)
    
    specify_i_range_manually = false
    if specify_i_range_manually
        i_range = [0] 
        # Manually specify the indices of the checkpoint files to be read in and processed e.g. i_range = [0, 1, 2, 3].
    else
        i_range = first_index:interval_index:last_index
    end

    for i in i_range
    
        filename = path * @sprintf("/neverworld_high_resolution_checkpoint_iteration%d.jld2", i)
        @printf("Extracting data from checkpoint file %s:\n", filename)
        T_Field, v_Field = CreateFieldsFromCheckPointOutput(grid, interior_size, halo_size, filename)
        
        int_T_xyz = compute!(Field(Integral(T_Field)))
        @printf("The heat content over the entire domain is %.6g.\n", int_T_xyz[1,1,1])
        iTimeSeries += 1
        int_T_xyz_TimeSeries[iTimeSeries] = int_T_xyz[1,1,1]
        
        v̄ = compute!(Field(Integral(v_Field, dims = 1)))
        ψ = Field{Nothing, Face, Face}(grid)
        
        for k in 2:grid.Nz+1
            dz = grid.Δzᵃᵃᶜ[k-1]      
            for j in 1:grid.Ny+1
                ψ[1, j, k] = ψ[1, j, k - 1] + v̄[1, j, k - 1] * dz
            end
        end
        
        ψ_Plot[:, :] = ψ[1, 1:grid.Ny+1, 1:grid.Nz+1]
        ψ_Mean_Plot[:, :] += ψ_Plot[:, :]    
        
        if make_animation
            ψ_Plot_TimeSeries[iTimeSeries, :, :] = ψ[1, 1:grid.Ny+1, 1:grid.Nz+1]
            title_TimeSeries[iTimeSeries] = @sprintf("Streamfunction Along Zonal Section at Output Index %d", i)        
        end
        
        if i == first_index || i == last_index
        
            if i == first_index
                title_Plot = "Initial Streamfunction Along Zonal Section"
                filename_Plot = "InitialStreamfunctionAlongZonalSection.pdf"
            elseif i == last_index
                title_Plot = "Final Streamfunction Along Zonal Section"
                filename_Plot = "FinalStreamfunctionAlongZonalSection.pdf"
            end
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                                     ψ_Plot[:, :], resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], 
                                     [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, "Streamfunction", 
                                     22.5, 10, 17.5, filename_Plot)
        
        end
        
    end
    
    ψ_Mean_Plot[:, :] /= n_indices
    MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                             ψ_Mean_Plot[:, :], resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], [17.5, 17.5], 
                             [10, 10], 1, "Time-Averaged Streamfunction Along Zonal Section", 27.5, 15, :balance, 100, 
                             "Streamfunction", 22.5, 10, 17.5, "TimeAveragedStreamfunctionAlongZonalSection.pdf")  
    
    if make_animation
        filename_Plot_Animation = "TimeEvolutionOfStreamfunctionAlongZonalSection"
        MakeHeatMapOrContourPlotAnimation(
        path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, ψ_Plot_TimeSeries[:, :, :], 
        resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], [17.5, 17.5], [10, 10], 1, title_TimeSeries, 27.5, 
        15, :balance, 100, "Streamfunction", 22.5, 10, 17.5, filename_Plot_Animation)
    end
    
    WriteOutputToFile1D(path, i_range, int_T_xyz_TimeSeries, "TimeEvolutionOfIntegratedHeatContent")
    indices, int_T_xyz_TimeSeries = ReadOutputFromFile1D(path, "TimeEvolutionOfIntegratedHeatContent.curve")
    MakeSingleLineOrScatterPlot(path, "scatter_line_plot", indices, int_T_xyz_TimeSeries, resolution, 2, :black, :rect,
                                0, ["Output Time Index", "Integrated Heat Content"], [25, 25], [17.5, 17.5], [10, 10], 
                                1, "Time Evolution of Integrated Heat Content", 27.5, 15, 
                                "TimeEvolutionOfIntegratedHeatContent.pdf")
        
end


function ComputeStreamFunctionAndPlotMeridionalOverturningCirculation_2(
path, specified_first_time_index, specified_last_time_index, use_all_time_indices = true)

    arch = CPU()
    new_degree = 1
    grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))
    
    filename_TData = path * "/Tdata.jld2"
    file_TData = jldopen(filename_TData)
    T_Array = file_TData["Tdata"]
    T_Field = Field{Center, Center, Center}(grid)
    
    filename_vData = path * "/Vdata.jld2"
    file_vData = jldopen(filename_vData)
    v_Array = file_vData["Vdata"]
    v_Field = Field{Center, Face, Center}(grid)
    
    interior_size = (grid.Nx, grid.Ny, grid.Nz) 
    halo_size = (grid.Hx, grid.Hy, grid.Hz)
    
    φᵃᶠᵃ_Array_Interior = ExtractInterior1DArray(grid.φᵃᶠᵃ, interior_size, halo_size, "afa")
    φᵃᶠᵃ_Array_Interior_Plot = zeros(grid.Ny+1)
    φᵃᶠᵃ_Array_Interior_Plot[:] = φᵃᶠᵃ_Array_Interior[:]
    
    zᵃᵃᶠ_Array_Interior = ExtractInterior1DArray(grid.zᵃᵃᶠ, interior_size, halo_size, "aaf")
    zᵃᵃᶠ_Array_Interior_Plot = zeros(grid.Nz+1)
    zᵃᵃᶠ_Array_Interior_Plot[:] = zᵃᵃᶠ_Array_Interior[:]/1000
    
    if use_all_time_indices == true
        first_time_index = 2 
        # Start from the second index since a zero velocity initial condition will result in a zero streamfunction, 
        # which in turn will throw an error when plotting the heat map or contour plot of the streamfunction.
        last_time_index = size(T_Array)[4]
    else
        first_time_index = specified_first_time_index
        last_time_index = specified_last_time_index
    end
    
    ψ_Plot = zeros(grid.Ny+1, grid.Nz+1)
    n_time_indices = last_time_index - first_time_index + 1
    ψ_Mean_Plot = zeros(grid.Ny+1, grid.Nz+1)
    
    make_animation = false
    if make_animation
        ψ_Plot_TimeSeries = zeros(n_time_indices, grid.Ny+1, grid.Nz+1)
        title_TimeSeries = fill("", n_time_indices)
    end
    
    int_T_xyz_TimeSeries = zeros(n_time_indices)
    iTimeSeries = 0
    
    resolution = (900, 750)
    
    for i in first_time_index:last_time_index

        @printf("Extracting data at time index %3d:\n", i)
        
        T_Array_Interior = ExtractInterior3DArray(T_Array[:, :, :, i], interior_size, halo_size, "ccc")
        v_Array_Interior = ExtractInterior3DArray(v_Array[:, :, :, i], interior_size, halo_size, "cfc")
        
        set!(T_Field, T_Array_Interior)
        set!(v_Field, v_Array_Interior)
        
        int_T_xyz = compute!(Field(Integral(T_Field)))
        @printf("The heat content over the entire domain is %.6g.\n", int_T_xyz[1,1,1])
        iTimeSeries += 1
        int_T_xyz_TimeSeries[iTimeSeries] = int_T_xyz[1,1,1]
        
        v̄ = compute!(Field(Integral(v_Field, dims = 1)))
        ψ = Field{Nothing, Face, Face}(grid)
        
        # Note that size(ψ) is (1, 141, 70) and size(ψ.data) is (1, 151, 80). 
        # This is because ψ.data is a 1×151×80 OffsetArray(::Array{Float64, 3}, 1:1, -4:146, -4:75).
        
        for k in 2:grid.Nz+1
            dz = grid.Δzᵃᵃᶜ[k-1]      
            for j in 1:grid.Ny+1
                ψ[1, j, k] = ψ[1, j, k - 1] + v̄[1, j, k - 1] * dz
            end
        end

        ψ_Plot[:, :] = ψ[1, 1:grid.Ny+1, 1:grid.Nz+1]
        ψ_Mean_Plot[:, :] += ψ_Plot[:, :]   
        
        if make_animation
            ψ_Plot_TimeSeries[iTimeSeries, :, :] = ψ[1, 1:grid.Ny+1, 1:grid.Nz+1]
            title_TimeSeries[iTimeSeries] = @sprintf("Streamfunction Along Zonal Section at Time Index %d", i)
        end 
        
        if i == first_time_index || i == last_time_index
        
            if i == first_time_index
                title_Plot = "Initial Streamfunction Along Zonal Section"
                filename_Plot = "InitialStreamfunctionAlongZonalSection.pdf"
            elseif i == last_time_index
                title_Plot = "Final Streamfunction Along Zonal Section"
                filename_Plot = "FinalStreamfunctionAlongZonalSection.pdf"
            end
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                                     ψ_Plot[:, :], resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], 
                                     [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, "Streamfunction", 
                                     22.5, 10, 17.5, filename_Plot)
            
        end
        
    end
    
    ψ_Mean_Plot[:, :] /= n_time_indices
    MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                             ψ_Mean_Plot[:, :], resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], [17.5, 17.5], 
                             [10, 10], 1, "Time-Averaged Streamfunction Along Zonal Section", 27.5, 15, :balance, 100, 
                             "Streamfunction", 22.5, 10, 17.5, "TimeAveragedStreamfunctionAlongZonalSection.pdf") 
    
    if make_animation
        filename_Plot_Animation = "TimeEvolutionOfStreamfunctionAlongZonalSection"
        MakeHeatMapOrContourPlotAnimation(
        path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, ψ_Plot_TimeSeries[:, :, :], 
        resolution, ["Latitude (degree)", "Depth (km)"], [25, 25], [17.5, 17.5], [10, 10], 1, title_TimeSeries, 27.5, 
        15, :balance, 100, "Streamfunction", 22.5, 10, 17.5, filename_Plot_Animation)
    end

    WriteOutputToFile1D(path, first_time_index:last_time_index, int_T_xyz_TimeSeries, 
                        "TimeEvolutionOfIntegratedHeatContent")
    time_indices, int_T_xyz_TimeSeries = ReadOutputFromFile1D(path, "TimeEvolutionOfIntegratedHeatContent.curve")
    MakeSingleLineOrScatterPlot(path, "scatter_line_plot", time_indices, int_T_xyz_TimeSeries, resolution, 2, :black, 
                                :rect, 0, ["Output Time Index", "Integrated Heat Content"], [25, 25], [17.5, 17.5], 
                                [10, 10], 1, "Time Evolution of Integrated Heat Content", 27.5, 15, 
                                "TimeEvolutionOfIntegratedHeatContent.pdf")
    
end

post_process_on_Satori = false
if post_process_on_Satori
    path = "/nobackup/users/sbishnu/WenoNeverworld_uq_of_bc_Output_Data"
else
    path = "../output"
end

Option = 2 # Choose Option to be 1 or 2. Default is 2.

if Option == 1

    first_index = 7200
    # Start from the second index since a zero velocity initial condition will result in a zero streamfunction, which in 
    # turn will throw an error when plotting the heat map or contour plot of the streamfunction.
    last_index = 28800 # On Satori, change the last index to 5256000.
    interval_index = 7200
    ComputeStreamFunctionAndPlotMeridionalOverturningCirculation_1(path, first_index, last_index, interval_index)
    
elseif Option == 2

    first_time_index = 2 
    # Start from the second index since a zero velocity initial condition will result in a zero streamfunction, which in 
    # turn will throw an error when plotting the heat map or contour plot of the streamfunction.
    last_time_index = 192
    use_all_time_indices = true
    ComputeStreamFunctionAndPlotMeridionalOverturningCirculation_2(path, first_time_index, last_time_index, 
                                                                   use_all_time_indices)
    
end