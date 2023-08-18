using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute
using GLMakie


dir = "/storage2/WenoNeverworldData/"
prefix    = "weno_eight_check"
variables = ("u", "v", "w", "b")
# stride = 3

# Load FieldTimeSeries of "u", "v", "w", "b" and compute integrated quantities
@info "loading files and calculating integrated variables"
neverworld_fields = all_fieldtimeseries(prefix, dir; checkpointer = true, variables, number_files = 20);
# kinetic_energy = Diagnostics.integral_kinetic_energy(neverworld_fields[:u], neverworld_fields[:b]; stride)
# integral_heat_content = Diagnostics.heat_content(neverworld_fields[:b]; stride)

grid  = neverworld_fields[:u].grid
times = neverworld_fields[:u].times
Nt = length(times)
iterations    = Nt-20:Nt
average_times = times[iterations]

u = neverworld_fields[:u]
v = neverworld_fields[:v]
w = neverworld_fields[:w]
b = neverworld_fields[:b]

# Average relevant variables
@info "Averaging relevant variables over" average_times
ū = Diagnostics.time_average(u, iterations)
v̄ = Diagnostics.time_average(v, iterations)
w_avg = Diagnostics.time_average(w, iterations)
b̄ = Diagnostics.time_average(b, iterations)
ψ = Diagnostics.calculate_eulerian_MOC(v̄)

@compute b̄ᵢ = Field(Average(b̄, dims = 1))

loc = location(u^2 + v^2)
EKE = FieldTimeSeries{loc...}(grid, average_times)
IKE = FieldTimeSeries{loc...}(grid, average_times)

@compute MKE = Field(0.5 * (ū^2 + v̄^2))

u′ = FieldTimeSeries{location(u)...}(grid, average_times)
v′ = FieldTimeSeries{location(v)...}(grid, average_times)
w′ = FieldTimeSeries{location(w)...}(grid, average_times)
b′ = FieldTimeSeries{location(b)...}(grid, average_times)

@info "Computing fluctuating fields" 
for (idx, iter) in enumerate(iterations)
    @info "at iter $iter of $(iterations[end])"
    set!(u′[idx], u[iter] - ū)
    set!(v′[idx], v[iter] - v̄)
    set!(w′[idx], w[iter] - w_avg)
    set!(b′[idx], b[iter] - b̄)

    set!(EKE[idx], 0.5 * (u′[idx]^2 + v′[idx]^2))
    set!(IKE[idx], 0.5 * (u[idx]^2 + v[idx]^2))
end

@compute MKE_avg = Field(Average(MKE, dims = (1, 3)))
@compute EKE_avg = Field(Average(EKE, dims = (1, 3)))

v′b′ = FieldTimeSeries{location(v[1] * b[1])...}(grid, average_times)
w′b′ = FieldTimeSeries{location(w[1] * b[1])...}(grid, average_times)

@info "Computing eddy fluxes" 
for (idx, iter) in enumerate(iterations)
    @info "at iter $iter of $(iterations[end])"
    set!(v′b′[idx], v′[iter] * b′[iter])
    set!(w′b′[idx], v′[iter] * b′[iter])
end

@info "Computation of residual circulation"
vb_avg = Diagnostics.time_average(v′b′)
wb_avg = Diagnostics.time_average(w′b′)

@compute vbᵢ_avg = Field(Average(vb_avg, dims = 1))
@compute wbᵢ_avg = Field(Average(wb_avg, dims = 1))

@compute ∂ybᵢ = Field(∂y(b̄ᵢ)) 
@compute ∂zbᵢ = Field(∂z(b̄ᵢ)) 

@compute μ = Field(- wbᵢ_avg / vbᵢ_avg * ∂zbᵢ/ ∂ybᵢ)


function plot_and_save_heatmap(data, filename, title, xlabel, ylabel; colorrange = (0, 0.02), level =10)
    fig = Figure(resolution=(2000, 1000))
    ax = Axis(fig[1, 1], xlabel=xlabel, xlabelsize=30, yticks=0:50:1200, xticklabelsize=30, ylabel=ylabel, ylabelsize=30, xticks=0:50:1120, yticklabelsize=30, title=title, titlesize=50)
    hm = GLMakie.heatmap!(ax, data, colorrange=colorrange, contours=true, levels = level)
    #GLMakie.contour!(ax, data, color=:black, linewidth=3, levels=level)
    display(fig)
    save("plotting/$filename.png", fig)
end

print("function defined succesfully")
# Plot and save heatmaps for each quantity
#plot_and_save_heatmap(interior(vbᵢ_avg,1, :, :), "time_averaged_eighth_moc", "Time-Averaged MOC 1/8∘", "latitude [∘]", "depth [m]", colorrange = (0, 0.02))

plot_and_save_heatmap(interior(IKE[20], :, :, 69), "time_averaged_eighth_IKE", "Time-Averaged IKE 1/8∘", "longitude [∘]", "latitude [∘]", colorrange = (0, 0.04))
plot_and_save_heatmap(interior(EKE[20], :, :, 69), "time_averaged_eighth_EKE", "Time-Averaged EKE 1/8∘", "longitude [∘]", "latitude [∘]", colorrange = (0, 0.04))
plot_and_save_heatmap(interior(MKE, :, :, 69), "time_averaged_MKE_eighth", "Mean Kinetic Energy 1/8∘", "longitude [∘]", "latitude [∘]", colorrange = (0, 0.04))


#plot_and_save_heatmap(interior(b, :, :, 69), "time_averaged_eighth_stratification", "Stratification 1/8∘", "latitude [∘]", "depth [m]", colorrange = (-0.04, 0.04))




#plot_and_save_heatmap(interior(u′[20], :, :, 69), "time_averaged_u", "Time-Averaged u", "latitude [∘]", "depth [m]", colorrange = (0, 0.02))
#plot_and_save_heatmap(interior(v′b′, :, :, 69), "v_prime_times_b_prime", "v' times b'", "latitude [∘]", "depth [m]", colorrange = (0, 0.02))
#plot_and_save_heatmap(interior(w′b′, :, :, 69), "w_prime_times_b_prime", "w' times b'", "latitude [∘]", "depth [m]", colorrange = (0, 0.02))
#plot_and_save_heatmap(interior(w_avg, :, :, 69), "time_averaged_w", "Time-Averaged w", "latitude [∘]", "depth [m]", colorrange = (0, 0.02))

using GLMakie
using JLD2, Oceananigans, Statistics      

# Load the data for the half-degree resolution
@info "Loading data..."
local_path = pwd()

hfile_1 = jldopen("/storage2/WenoNeverworldData/weno_eight_checkpoint_iteration22782.jld2 ", "r")
oceangrid_1 = hfile_1["grid"]
halo = 7
z = oceangrid_1.underlying_grid.zᵃᵃᶜ[1:end-halo]
lat = collect(oceangrid_1.underlying_grid.φᵃᶜᵃ[1:end-halo])
lon_index = round(Int, size(b)[1]/2)
fig = Figure(resolution=(2000, 1000))
ax = Axis(fig[1, 1], xlabel="Latitude [∘]", xlabelsize=30, yticks=-4000:1000:0, xticklabelsize=30, ylabel="Depth [m]", ylabelsize=30, xticks=-90:20:1120, yticklabelsize=30, title="Heatmap of Time-Averaged b 1/8", titlesize=50)
hm = GLMakie.heatmap!(ax, lat, z, interior(b, lon_index, :, :), colormap= :plasma, contours = true, levels =10)  #colorrange = (7.5e-5, 3.6e-5))
contour!(ax, lat, z, interior(b, lon_index, :, :), color=:black, linewidth=3, levels=10, labels = true,
labelsize = 30, labelfont = :bold, labelcolor = :black)
display(fig)
save("plotting/moc_eight_final.png", fig)
