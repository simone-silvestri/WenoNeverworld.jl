using FilePaths
using CairoMakie


function MakeSingleLineOrScatterPlot(output_directory, plot_type, x, y, resolution, linewidth, linecolor, marker, 
                                     markersize, labels, labelsizes, ticklabelsizes, labelpaddings, aspect, title, 
                                     titlesize, titlegap, file_name)
    
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path) 
        mkdir(path) 
    end
    cd(path)
    
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1,1]; xlabel = labels[1], ylabel = labels[2], xlabelsize = labelsizes[1], ylabelsize = labelsizes[2], 
              xticklabelsize = ticklabelsizes[1], yticklabelsize = ticklabelsizes[2], xlabelpadding = labelpaddings[1], 
              ylabelpadding = labelpaddings[2], aspect = aspect, title = title, titlesize = titlesize, 
              titlegap = titlegap, titlefont = :bold)
    
    if plot_type == "line_plot"
        lines!(ax, x, y, linewidth = linewidth, color = linecolor)
    elseif plot_type == "scatter_plot"
        scatter!(ax, x, y, marker = marker, markersize = markersize, color = linecolor)
    elseif plot_type == "scatter_line_plot"
        scatterlines!(ax, x, y, linewidth = linewidth, marker = marker, markersize = markersize, color = linecolor)
    end
    
    save(file_name, fig)
    cd(cwd)
    
end


function MakeHeatMapOrContourPlot(output_directory, plot_type, x, y, φ, resolution, labels, labelsizes, ticklabelsizes, 
                                  labelpaddings, aspect, title, titlesize, titlegap, φ_limits, colormap, contour_levels, 
                                  file_name, specify_axis_limits = true, use_specified_limits = false)
    
    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path) 
        mkdir(path) 
    end
    cd(path)
    
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1,1]; xlabel = labels[1], ylabel = labels[2], xlabelsize = labelsizes[1], ylabelsize = labelsizes[2], 
    xticklabelsize = ticklabelsizes[1], yticklabelsize = ticklabelsizes[2], xlabelpadding = labelpaddings[1], 
    ylabelpadding = labelpaddings[2], aspect = aspect, title = title, titlesize = titlesize, 
    titlegap = titlegap, titlefont = :bold)
    
    if specify_axis_limits
        xlims!(ax, (minimum(x), maximum(x)))
        ylims!(ax, (minimum(y), maximum(y)))
    end
    
    if !use_specified_limits
        φ_limits = [minimum(φ), maximum(φ)]
    end
    
    if plot_type == "heat_map"
        hm = heatmap!(ax, x, y, φ; colorrange = φ_limits, colormap = colormap)
    elseif plot_type == "filled_contour_plot"
        hm = contourf!(ax, x, y, φ; levels = range(φ_limits..., length=contour_levels), colormap = colormap)  
    end
    
    Colorbar(fig[1,2], hm)
    
    save(file_name, fig)
    cd(cwd)
    
end


TestMakePlots = false

if TestMakePlots

    x = range(0, stop = 2π, length = 100)
    y = sin.(x)
    resolution = (800, 750)
    
    file_name = "LinePlotExample.pdf"
    MakeSingleLineOrScatterPlot("../output", "line_plot", x, y, resolution, 2, :black, :rect, 10, ["x", "sin(x)"], 
                                [25, 25], [17.5, 17.5], [10, 10], 1, "sin(x) vs x", 27.5, 15, file_name)
    
    file_name = "ScatterPlotExample.pdf"
    MakeSingleLineOrScatterPlot("../output", "scatter_plot", x, y, resolution, 2, :black, :rect, 10, ["x", "sin(x)"], 
                                [25, 25], [17.5, 17.5], [10, 10], 1, "sin(x) vs x", 27.5, 15, file_name)
    
    file_name = "ScatterLinePlotExample.pdf"
    MakeSingleLineOrScatterPlot("../output", "scatter_line_plot", x, y, resolution, 2, :black, :rect, 10, 
                                ["x", "sin(x)"], [25, 25], [17.5, 17.5], [10, 10], 1, "sin(x) vs x", 27.5, 15, 
                                file_name)
    
    nX = 100
    nY = 100
    x = range(0, stop = 2π, length = nX)
    y = range(0, stop = 2π, length = nY)
    φ = fill(0.0, nX, nY)
    for i in 1:nX
        for j in 1:nY
            φ[i,j] = sin(x[i]) * sin(y[j])
        end
    end
    resolution = (800, 750)
    φ_limits = [-1, 1]
    
    file_name = "HeatMapExample.pdf"
    MakeHeatMapOrContourPlot("../output", "heat_map", x, y, φ, resolution, ["x", "y"], [25, 25], [17.5, 17.5], 
                             [10, 10], 1, "sin(x) vs x", 27.5, 15, φ_limits, :balance, 100, file_name)
    
    file_name = "FilledContourPlotExample.pdf"
    MakeHeatMapOrContourPlot("../output", "filled_contour_plot", x, y, φ, resolution, ["x", "y"], [25, 25], 
                             [17.5, 17.5], [10, 10], 1, "sin(x) vs x", 27.5, 15, φ_limits, :balance, 100, file_name)
    
end