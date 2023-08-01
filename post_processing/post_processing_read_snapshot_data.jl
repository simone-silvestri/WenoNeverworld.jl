using Pkg
using Oceananigans
using JLD2
using CairoMakie
using Printf


include("post_processing_make_plots.jl")


function ReadSnapshotDataAndPlotIntegratedHeatContent()

    file_name = "../output/uq_one_degree_salinity_snapshots.jld2" 
    file = jldopen(file_name)
    iterations = parse.(Int,keys(file["timeseries/t"]))
    last_index = length(iterations)
    
    T_series = FieldTimeSeries(file_name, "T"; architecture = CPU())

    xT, yT, zT = nodes(T_series[1])
    
    T = T_series[last_index]
    
    # Compute the integral of temperature in the vertical direction.
    
    int_T_z = compute!(Field(Integral(T, dims = 3)))
    int_T_z_interior = interior(int_T_z, :, :, size(int_T_z, 3)) 
    
    # Compute the integral of temperature in all directions.
    
    int_T_xyz = compute!(Field(Integral(T)))
    
    @printf("The heat content over the entire domain is %.6g.", int_T_xyz[1,1,1])
    
    # Plot the depth integrated heat content.
    MakeHeatMapOrContourPlot("../output", "filled_contour_plot", xT, yT, int_T_z_interior, (850, 750), 
                             ["Latitude (degree)", "Longitude (degree)"], [25, 25], [17.5, 17.5], [10, 10], 1, 
                             "Depth Integrated Heat Content", 27.5, 15, :balance, 100, "Heat content", 22.5, 10, 17.5,
                             "DepthIntegratedHeatContent.pdf"; make_colorbar_symmetric_about_zero = true)
    
    nothing # hide
    
end


ReadSnapshotDataAndPlotIntegratedHeatContent()