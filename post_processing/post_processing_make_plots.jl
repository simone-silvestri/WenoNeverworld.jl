using FilePaths
using DelimitedFiles
using CairoMakie
using Printf


function WriteOutputToFile1D(output_directory, x, y, filename)

    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)
    
    filename *= ".curve"
    outputfile = open(filename, "w")
    
    write(outputfile, "#phi\n")
    for i in eachindex(x)
        write(outputfile, string(x[i], " ", y[i], "\n"))
        # The above line is equivalent to println(outputfile, "$(x[i]) $(y[i])")
    end
    
    close(outputfile)
    cd(cwd)

end


function ReadOutputFromFile1D(output_directory, filename)

    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path)
        mkdir(path)
    end
    cd(path)
    
    data = []
    count = 1
    open(filename, "r") do infile
        for line in eachline(infile)
            if count != 1
                push!(data, line)
            end
            count += 1
        end
    end
    data = readdlm(IOBuffer(join(data, "\n")))
    
    N = size(data, 1)
    x = zeros(N)
    y = zeros(N)
    for i in 1:N
        x[i] = data[i,1]
        y[i] = data[i,2]
    end
    
    cd(cwd)
    
    return (x, y)
    
end


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
                                  labelpaddings, aspect, title, titlesize, titlegap, colormap, contour_levels, cb_label, 
                                  cb_labelsize, cb_labelpadding, cb_ticksize, file_name, specify_axis_limits = true, 
                                  use_specified_φ_limits = false, specified_φ_limits = [0, 0], 
                                  make_colorbar_symmetric_about_zero = false, extrema_reduction_factor = 1)
    
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
    
    if use_specified_φ_limits
        φ_limits = specified_φ_limits
    else
        if make_colorbar_symmetric_about_zero
            φ_limits = [-maximum(abs.(φ)), maximum(abs.(φ))]
        else
            φ_limits = [minimum(φ), maximum(φ)]
        end
    end  
    φ_limits *= extrema_reduction_factor
    
    if plot_type == "heat_map"
        hm = heatmap!(ax, x, y, φ; colorrange = φ_limits, colormap = colormap)
    elseif plot_type == "filled_contour_plot"
        hm = contourf!(ax, x, y, φ; levels = range(φ_limits..., length=contour_levels), colormap = colormap)  
    end
    
    Colorbar(fig[1,2], hm; label = cb_label, labelsize = cb_labelsize, labelpadding = cb_labelpadding, 
             ticksize = cb_ticksize)
    
    save(file_name, fig)
    cd(cwd)
    
end


function MakeHeatMapOrContourPlotAnimation(
output_directory, plot_type, x, y, φ_time_series, resolution, labels, labelsizes, ticklabelsizes, labelpaddings, aspect, 
title_time_series, titlesize, titlegap, colormap, contour_levels, cb_label, cb_labelsize, cb_labelpadding, 
cb_ticksize, file_name, specify_axis_limits = true, use_specified_φ_limits = false, specified_φ_limits = [0, 0], 
make_colorbar_symmetric_about_zero = false, extrema_reduction_factor = 1, frame_rate = 10)

    cwd = pwd()
    path = joinpath(cwd, output_directory)
    if !isdir(path) 
        mkdir(path) 
    end
    cd(path)

    if use_specified_φ_limits
        φ_limits = specified_φ_limits
    else
        if make_colorbar_symmetric_about_zero
            φ_limits = [-maximum(abs.(φ_time_series)), maximum(abs.(φ_time_series))]
        else
            φ_limits = [minimum(φ_time_series), maximum(φ_time_series)]
        end
    end  
    φ_limits *= extrema_reduction_factor
    
    @info "Creating animation"
    
    n = Observable(1)
    φ = @lift φ_time_series[$n, :, :]
    fig = Figure(resolution = resolution)
    title = @lift title_time_series[$n]
    ax = Axis(fig[1,1]; xlabel = labels[1], ylabel = labels[2], xlabelsize = labelsizes[1], ylabelsize = labelsizes[2], 
              xticklabelsize = ticklabelsizes[1], yticklabelsize = ticklabelsizes[2], xlabelpadding = labelpaddings[1], 
              ylabelpadding = labelpaddings[2], aspect = aspect, title = title, titlesize = titlesize, 
              titlegap = titlegap, titlefont = :bold)    
    if specify_axis_limits
        xlims!(ax, (minimum(x), maximum(x)))
        ylims!(ax, (minimum(y), maximum(y)))
    end
    
    if plot_type == "heat_map"
        hm = heatmap!(ax, x, y, φ; colorrange = φ_limits, colormap = colormap)
    elseif plot_type == "filled_contour_plot"
        hm = contourf!(ax, x, y, φ; levels = range(φ_limits..., length=contour_levels), colormap = colormap)  
    end
    
    Colorbar(fig[1,2], hm; label = cb_label, labelsize = cb_labelsize, labelpadding = cb_labelpadding, 
             ticksize = cb_ticksize)
    
    frames = 1:size(φ_time_series)[1]
    
    CairoMakie.record(fig, file_name * ".mp4", frames, framerate = frame_rate) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[] = i
    end
    
    cd(cwd)
    
    nothing # hide output

end


TestWriteOutputToAndReadOutputFromFile1D = false

if TestWriteOutputToAndReadOutputFromFile1D

    x = collect(0.0:1.0:9.0)
    y = collect(0.0:2.0:18.0)
    output_directory = "../output/"
    WriteOutputToFile1D(output_directory, x, y, "TestWriteOutputToFile1D")
    
    x1 = collect(0.0:1.0:9.0)
    y1 = collect(0.0:2.0:18.0)
    println("Write to file:")
    for i in eachindex(x1)
        println(x1[i], " ", y1[i])
    end
    
    output_directory = "../output/"
    WriteOutputToFile1D(output_directory, x1, y1, "TestWriteOutputToFile1D")
    x2, y2 = ReadOutputFromFile1D(output_directory, "TestWriteOutputToFile1D.curve")
    
    println(" ")
    println("Read from file:")
    for i in eachindex(x2)
        println(x2[i], " ", y2[i])
    end
    
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
    nTime = 100
    t = range(0, stop = 2π, length = nTime)
    φ_time_series = fill(0.0, nTime, nX, nY)
    title_time_series = fill("", nTime)
    
    for iX in 1:nX
        for iY in 1:nY
            φ[iX, iY] = sin(x[iX]) * sin(y[iY])
            for iTime in 1:nTime
                φ_time_series[iTime, iX, iY] = sin(x[iX]) * sin(y[iY]) * sin(t[iTime])
                title_time_series[iTime] = @sprintf("sin(x) * sin(y) * sin(t) at t = %.2f", t[iTime])
            end
        end
    end
    
    resolution = (850, 750)
    
    file_name = "HeatMapExample.pdf"
    MakeHeatMapOrContourPlot("../output", "heat_map", x, y, φ, resolution, ["x", "y"], [25, 25], [17.5, 17.5], 
                             [10, 10], 1, "sin(x) * sin(y)", 27.5, 15, :balance, 100, "Heat map colorbar", 22.5, 
                             10, 17.5, file_name)
    
    file_name = "FilledContourPlotExample.pdf"
    MakeHeatMapOrContourPlot("../output", "filled_contour_plot", x, y, φ, resolution, ["x", "y"], [25, 25], 
                             [17.5, 17.5], [10, 10], 1, "sin(x) * sin(y)", 27.5, 15, :balance, 100, 
                             "Filled contour plot colorbar", 22.5, 10, 17.5, file_name)
    
    file_name = "HeatMapExampleAnimation"
    MakeHeatMapOrContourPlotAnimation(
    "../output", "heat_map", x, y, φ_time_series, resolution, ["x", "y"], [25, 25], [17.5, 17.5], [10, 10], 1, 
    title_time_series, 27.5, 15, :balance, 100, "Heat map colorbar", 22.5, 10, 17.5, file_name)
    
    file_name = "FilledContourPlotExampleAnimation"
    MakeHeatMapOrContourPlotAnimation(
    "../output", "filled_contour_plot", x, y, φ_time_series, resolution, ["x", "y"], [25, 25], [17.5, 17.5], [10, 10], 
    1, title_time_series, 27.5, 15, :balance, 100, "Filled contour plot colorbar", 22.5, 10, 17.5, file_name)
    
end