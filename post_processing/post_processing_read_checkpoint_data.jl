using Pkg
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.AbstractOperations: KernelFunctionOperation
using JLD2
using CairoMakie
using Printf
using WenoNeverworld


include("post_processing_make_plots.jl")


function ComputeStreamFunctionAndPlotMeridionalOverturningCirculation(path, first_index, last_index, range_of_indices, 
                                                                      n_indices)

    arch = CPU()
    new_degree = 0.25
    grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70))
    
    Nx, Ny, Nz = grid.Nx, grid.Ny, grid.Nz
    Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz
    
    λᶜᵃᵃ_Array_Interior_Plot = grid.λᶜᵃᵃ[1:Nx]
    λᶠᵃᵃ_Array_Interior_Plot_Periodic = grid.λᶠᵃᵃ[1:Nx]
    φᵃᶜᵃ_Array_Interior_Plot = grid.φᵃᶜᵃ[1:Ny]
    φᵃᶠᵃ_Array_Interior_Plot = grid.φᵃᶠᵃ[1:Ny+1]
    zᵃᵃᶜ_Array_Interior_Plot = grid.zᵃᵃᶜ[1:Nz]/1000
    zᵃᵃᶠ_Array_Interior_Plot = grid.zᵃᵃᶠ[1:Nz+1]/1000
    
    ψ_Plot_AlongZonalSection = zeros(Ny+1, Nz+1)
    ψ_Mean_Plot_AlongZonalSection = zeros(Ny+1, Nz+1)
    
    int_T_xyz_TimeSeries = zeros(n_indices)
    iTimeSeries = 0
    
    resolution = (900, 750)
    
    for i in range_of_indices
    
        filename = path * @sprintf("/neverworld_seawater_quarter_checkpoint_iteration%d.jld2", i)
        @printf("Extracting data from checkpoint file %s:\n", filename)
        
        file = jldopen(filename)
        
        T_Field = CenterField(grid)
        v_Field = YFaceField(grid)
        
        set!(T_Field, file["T/data"][1 + Hx : Nx + Hx, 1 + Hy : Ny + Hy, 1 + Hz : Nz + Hz])
        set!(v_Field, file["v/data"][1 + Hx : Nx + Hx, 1 + Hy : Ny + Hy + 1, 1 + Hz : Nz + Hz])
        
        #=
        The above lines are equivalent to the following lines:
        set!(T_Field, file["T/data"][1 + Hx : end - Hx, 1 + Hy : end - Hy, 1 + Hz : end - Hz])
        set!(v_Field, file["v/data"][1 + Hx : end - Hx, 1 + Hy : end - Hy, 1 + Hz : end - Hz])
        =#
        
        if i == last_index
        
            S_Field = CenterField(grid)
            u_Field = XFaceField(grid)
            
            set!(S_Field, file["S/data"][1 + Hx : Nx + Hx, 1 + Hy : Ny + Hy, 1 + Hz : Nz + Hz])
            set!(u_Field, file["u/data"][1 + Hx : Nx + Hx, 1 + Hy : Ny + Hy, 1 + Hz : Nz + Hz])
            
            #=
            The above lines are equivalent to the following lines:
            set!(S_Field, file["S/data"][1 + Hx : end - Hx, 1 + Hy : end - Hy, 1 + Hz : end - Hz])
            set!(u_Field, file["u/data"][1 + Hx : end - Hx, 1 + Hy : end - Hy, 1 + Hz : end - Hz])
            =#
            
        end
        
        int_T_xyz = compute!(Field(Integral(T_Field)))
        @printf("The heat content over the entire domain is %.6g.\n", int_T_xyz[1,1,1])
        iTimeSeries += 1
        int_T_xyz_TimeSeries[iTimeSeries] = int_T_xyz[1,1,1]

        v̄_AlongZonalSection = compute!(Field(Average(v_Field, dims = 1)))
        
        if i == last_index
            
            T_Plot_OnSurface = interior(T_Field, :, :, Nz)
            S_Plot_OnSurface = interior(S_Field, :, :, Nz)
            u_Plot_OnSurface = interior(u_Field, :, :, Nz)
            v_Plot_OnSurface = interior(v_Field, :, :, Nz)
            
            vorticity_operation = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u_Field, v_Field)
            ζ_Field = compute!(Field(vorticity_operation))
            ζ_Plot_OnSurface = interior(ζ_Field, :, :, Nz)
            
            T̄_AlongZonalSection = compute!(Field(Average(T_Field, dims = 1)))
            S̄_AlongZonalSection = compute!(Field(Average(S_Field, dims = 1)))
            ū_AlongZonalSection = compute!(Field(Average(u_Field, dims = 1)))
            
            T_Plot_AlongZonalSection = interior(T̄_AlongZonalSection, 1, :, :)
            S_Plot_AlongZonalSection = interior(S̄_AlongZonalSection, 1, :, :)
            u_Plot_AlongZonalSection = interior(ū_AlongZonalSection, 1, :, :)
            v_Plot_AlongZonalSection = interior(v̄_AlongZonalSection, 1, :, :)    
            
        end
        
        ψ_AlongZonalSection = Field{Nothing, Face, Face}(grid)
        
        for k in 2:Nz+1
            dz = grid.Δzᵃᵃᶜ[k-1]      
            for j in 1:Ny+1
                ψ_AlongZonalSection[1, j, k] = ψ_AlongZonalSection[1, j, k - 1] + v̄_AlongZonalSection[1, j, k - 1] * dz
            end
        end
        
        ψ_Plot_AlongZonalSection[:, :] = ψ_AlongZonalSection[1, 1:Ny+1, 1:Nz+1]
        ψ_Mean_Plot_AlongZonalSection[:, :] += ψ_Plot_AlongZonalSection[:, :]    
        
        if i == last_index
            
            # Plots on the surface
            
            title_Plot = "Final Temperature on the Surface"
            filename_Plot = "FinalTemperatureOnSurface.png"
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", λᶜᵃᵃ_Array_Interior_Plot,  φᵃᶜᵃ_Array_Interior_Plot, 
                                     T_Plot_OnSurface[:, :], resolution, ["Latitude (degree)", "Longitude (degree)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Temperature", 22.5, 10, 17.5, filename_Plot)
            
            title_Plot = "Final Salinity on the Surface"
            filename_Plot = "FinalSalinityOnSurface.png"
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", λᶜᵃᵃ_Array_Interior_Plot,  φᵃᶜᵃ_Array_Interior_Plot, 
                                     S_Plot_OnSurface[:, :], resolution, ["Latitude (degree)", "Longitude (degree)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Salinity", 22.5, 10, 17.5, filename_Plot)

            title_Plot = "Final Zonal Velocity on the Surface"
            filename_Plot = "FinalZonalVelocityOnSurface.png"
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", λᶠᵃᵃ_Array_Interior_Plot_Periodic, 
                                     φᵃᶜᵃ_Array_Interior_Plot, u_Plot_OnSurface[:, :], resolution, 
                                     ["Latitude (degree)", "Longitude (degree)"], [25, 25], [17.5, 17.5], [10, 10], 1, 
                                     title_Plot, 27.5, 15, :balance, 100, "Zonal Velocity", 22.5, 10, 17.5, 
                                     filename_Plot)

            title_Plot = "Final Meridional Velocity on the Surface"
            filename_Plot = "FinalMeridionalVelocityOnSurface.png"
                                    
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", λᶜᵃᵃ_Array_Interior_Plot,  φᵃᶠᵃ_Array_Interior_Plot, 
                                     v_Plot_OnSurface[:, :], resolution, ["Latitude (degree)", "Longitude (degree)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Meridional Velocity", 22.5, 10, 17.5, filename_Plot)
            
            title_Plot = "Final Relative Vorticity on the Surface"
            filename_Plot = "FinalRelativeVorticityOnSurface.png"
                                    
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", λᶠᵃᵃ_Array_Interior_Plot_Periodic,  
                                     φᵃᶠᵃ_Array_Interior_Plot, ζ_Plot_OnSurface[:, :], resolution, 
                                     ["Latitude (degree)", "Longitude (degree)"], [25, 25], [17.5, 17.5], [10, 10], 1, 
                                     title_Plot, 27.5, 15, :balance, 100, "Relative Vorticity", 22.5, 10, 17.5, 
                                     filename_Plot)            
            
            # Zonally averaged plots
            
            title_Plot = "Final Temperature Along Zonal Section"
            filename_Plot = "FinalTemperatureAlongZonalSection.png"
            n_levels = 5
            d_level = (maximum(T_Plot_AlongZonalSection[:, :]) - minimum(T_Plot_AlongZonalSection[:, :]))/n_levels
            levels = minimum(T_Plot_AlongZonalSection[:, :]) : d_level : maximum(T_Plot_AlongZonalSection[:, :])
            
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶜᵃ_Array_Interior_Plot, zᵃᵃᶜ_Array_Interior_Plot, 
                                     T_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Temperature", 22.5, 10, 17.5, filename_Plot, plot_contour_lines = true, 
                                     levels = levels)
            
            title_Plot = "Final Salinity Along Zonal Section"
            filename_Plot = "FinalSalinityAlongZonalSection.png"
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶜᵃ_Array_Interior_Plot, zᵃᵃᶜ_Array_Interior_Plot, 
                                     S_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Salinity", 22.5, 10, 17.5, filename_Plot)
            
            title_Plot = "Final Zonal Velocity Along Zonal Section"
            filename_Plot = "FinalZonalVelocityAlongZonalSection.png"
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶜᵃ_Array_Interior_Plot, zᵃᵃᶜ_Array_Interior_Plot, 
                                     u_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Zonal Velocity", 22.5, 10, 17.5, filename_Plot)
            
            title_Plot = "Final Meridional Velocity Along Zonal Section"
            filename_Plot = "FinalMeridionalVelocityAlongZonalSection.png"
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶜ_Array_Interior_Plot, 
                                     v_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Meridional Velocity", 22.5, 10, 17.5, filename_Plot)
            
            title_Plot = "Final Streamfunction Along Zonal Section"
            filename_Plot = "FinalStreamfunctionAlongZonalSection.png"
            n_levels = 5
            d_level = (maximum(ψ_Plot_AlongZonalSection[:, :]) - minimum(ψ_Plot_AlongZonalSection[:, :]))/n_levels
            levels = minimum(ψ_Plot_AlongZonalSection[:, :]) : d_level : maximum(ψ_Plot_AlongZonalSection[:, :])
            MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                                     ψ_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                                     [25, 25], [17.5, 17.5], [10, 10], 1, title_Plot, 27.5, 15, :balance, 100, 
                                     "Streamfunction", 22.5, 10, 17.5, filename_Plot; 
                                     make_colorbar_symmetric_about_zero = true, plot_contour_lines = true, 
                                     levels = levels)
            
        end

    end
    
    ψ_Mean_Plot_AlongZonalSection[:, :] /= n_indices
    n_levels = 5
    d_level = (maximum(ψ_Mean_Plot_AlongZonalSection[:, :]) - minimum(ψ_Mean_Plot_AlongZonalSection[:, :]))/n_levels
    levels = minimum(ψ_Mean_Plot_AlongZonalSection[:, :]) : d_level : maximum(ψ_Mean_Plot_AlongZonalSection[:, :])
    MakeHeatMapOrContourPlot(path, "filled_contour_plot", φᵃᶠᵃ_Array_Interior_Plot, zᵃᵃᶠ_Array_Interior_Plot, 
                             ψ_Mean_Plot_AlongZonalSection[:, :], resolution, ["Latitude (degree)", "Depth (km)"], 
                             [25, 25], [17.5, 17.5], [10, 10], 1, "Time-Averaged Streamfunction Along Zonal Section", 
                             27.5, 15, :balance, 100, "Streamfunction", 22.5, 10, 17.5, 
                             "TimeAveragedStreamfunctionAlongZonalSection.png"; 
                             make_colorbar_symmetric_about_zero = true, plot_contour_lines = true, levels = levels)
    
    WriteOutputToFile1D(path, range_of_indices, int_T_xyz_TimeSeries, "TimeEvolutionOfIntegratedHeatContent")
    indices, int_T_xyz_TimeSeries = ReadOutputFromFile1D(path, "TimeEvolutionOfIntegratedHeatContent.curve")
    MakeSingleLineOrScatterPlot(path, "scatter_line_plot", indices, int_T_xyz_TimeSeries, resolution, 2, :black, :rect,
                                0, ["Output Time Index", "Integrated Heat Content"], [25, 25], [17.5, 17.5], [10, 10], 
                                1, "Time Evolution of Integrated Heat Content", 27.5, 15, 
                                "TimeEvolutionOfIntegratedHeatContent.png")
        
end


post_process_on_Satori = true
if post_process_on_Satori
    path = "/nobackup/users/sbishnu/WenoNeverworld_uq_of_bc_Output_Data"
else
    path = "../output"
end


specify_range_of_indices_manually = true
if specify_range_of_indices_manually
    range_of_indices = [137373, 247202, 351634, 456529]
    # Manually specify the indices of the checkpoint files to be read in and processed 
    # e.g. range_of_indices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] for 10 checkpoint files. 
    #
    # Also, don't start from 0 since a zero velocity initial condition will result in a zero streamfunction, which in 
    # turn will throw an error when plotting the heat map or contour plot of the streamfunction.
    first_index = range_of_indices[1]
    last_index = range_of_indices[end]
else
    first_index = 7200
    # Start from the second index since a zero velocity initial condition will result in a zero streamfunction, which in
    # turn will throw an error when plotting the heat map or contour plot of the streamfunction.
    last_index = 28800 # On Satori, change the last index to 5256000.
    interval_index = 7200
    range_of_indices = first_index:interval_index:last_index
end
n_indices = length(range_of_indices)
ComputeStreamFunctionAndPlotMeridionalOverturningCirculation(path, first_index, last_index, range_of_indices, n_indices)