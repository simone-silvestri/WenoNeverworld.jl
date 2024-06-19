using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Fields: condition_operand
using Oceananigans.AbstractOperations: materialize_condition!
using GLMakie
using Statistics
using JLD2


f = all_fieldtimeseries("weno_thirtytwo_compressed_iteration_new1507554.jld2", "/storage4/WenoNeverworldData/"; checkpointer = true)

grid = f[:u].grid
b = f[:b][1];
Nx, Ny, Nz = size(grid)

u = f[:u][1];
v = f[:v][1];

KE = KineticEnergy(f, 1; indices = (:, :, 69));
KE_field = Field((Center, Center, Nothing), grid.underlying_grid);
set!(KE_field, interior(KE, :, :, 1));

ζ  = VerticalVorticity(f, 1; indices = (:, :, 69));
Ld = DeformationRadius(f, 1);
N² = Stratification(f, 1; indices = (:, :, 69));
N²_mid = Stratification(f, 1; indices = (:, :, 45));
q  = PotentialVorticity(f, 1; indices = (:, :, 69));

#=
f1 = all_fieldtimeseries("neverworld_backscatter_checkpoint_iteration288000.jld2"; checkpointer = true)
grid1 = f1[:u].grid
u1 = f1[:u][1];
v1 = f1[:v][1];
KE1 = KineticEnergy(f1, 1; indices = (:, :, 69));
KE1_field = Field((Center, Center, Nothing), grid1.underlying_grid);
set!(KE1_field, interior(KE1, :, :, 1));

KE_coarse = similar(KE1_field);
=#
#WenoNeverworld.Auxiliaries.three_dimensional_regrid!(KE_coarse, KE_field);
#WenoNeverworld.Auxiliaries.three_dimensional_regrid!(KE_field);

nothing

xF, yF, _ = nodes(ζ)
xC, yC, _ = nodes(KE)

xF = xF .- 30
xC = xC .- 30

#=
fig = Figure(size = (1000, 1000), fontsize = 15)
ax  = Axis(fig[1, 1],
          ylabel = L"\text{Latitude}",
          xlabel = L"\text{Longitude}",
          xticks = ([0, 15, 30, 45, 60], [L"0", L"15", L"30", L"45", L"60"]),
          yticks = ([-70, -35, 0, 35, 70], [L"-70", L"-35", L"0", L"30", L"70"]))

hm1 = heatmap!(ax, xC, yC, log10.(interior(KE, :, :, 1)), colormap = :magma, colorrange = (-3, 0))
cb  = Colorbar(fig[0, 1], hm1, vertical = false, label = L"\text{Surface Kinetic Energy [m}^2\text{s}^{-2}\text{]}",
               ticks = ([-3, -2, -1, 0], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{-0}"]))
hidedecorations!(ax)
hidespines!(ax)
display(fig)
=#

fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/32∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:20:70, yticksize = 15, xticksize = 15)
        
hm1 = heatmap!(ax, xF, yF, interior(ζ, :, :, 1), colormap = :berlin, colorrange = (-6e-5, 6e-5))
#cb  = Colorbar(fig[0, 2], hm1, vertical = false, label = L"\text{Surface Vertical Vorticity [s}^{-1} \cdot 10^{-5}\text{]}", ticks = ([-3e-5, 0, 3e-5], [L"-3", L"0", L"3"]))
cb = Colorbar(fig[1,2], hm1, width = 30, ticksize = 10, ticklabelsize = 40, ticks = ([-3e-5, 0, 3e-5], [L"-3", L"0", L"3"]), height = Relative(3/4))
display(fig)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("vort_32_cb.png", fig, px_per_unit = 5)

#=
fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/32∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:20:70, yticksize = 15, xticksize = 15)

hm1 = heatmap!(ax, xC, yC, interior(b, :, :, 69) ./ 2e-3, colormap = :thermal, colorrange = (0, 30))
#cb  = Colorbar(fig[0, 3], hm1, vertical = false, label = L"\text{Surface Temperature [}^\circ\text{C}^{-1}\text{]}", ticks = ([0, 10, 20, 30], [L"0", L"10", L"20", L"30"]))
cb = Colorbar(fig[1,2], hm1, width = 30, ticksize = 10, ticklabelsize = 40, height = Relative(3/4))
display(fig)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("surf_temp_32.png", fig, px_per_unit = 5)
=#