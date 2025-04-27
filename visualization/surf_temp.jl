using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Fields: condition_operand
using Oceananigans.AbstractOperations: materialize_condition!
using CairoMakie
using Statistics
using JLD2
using FFTW


f = all_fieldtimeseries("weno_eighth_checkpoint_iteration28667520.jld2", "/storage4/WenoNeverworldData/sixteenth_degree_new/"; checkpointer = true)

#weno_quarter__checkpoint_iteration70005600.jld2
#weno_sixteen_checkpoint_iteration1156320.jld2
#weno_eighth_checkpoint_iteration28667520.jld2
#weno_half_checkpoint_iteration42167520.jld2

grid = f[:u].grid
b = f[:b][1];
Nx, Ny, Nz = size(grid)

u = f[:u][1];
v = f[:v][1];

KE = KineticEnergy(f, 1; indices = (:, :, 34));
KE_field = Field((Center, Center, Nothing), grid.underlying_grid);
set!(KE_field, interior(KE, :, :, 1));

ζ  = VerticalVorticity(f, 1; indices = (:, :, 34));
Ld = DeformationRadius(f, 1);
N² = Stratification(f, 1; indices = (:, :, 34));
N²_mid = Stratification(f, 1; indices = (:, :, 25));
q  = PotentialVorticity(f, 1; indices = (:, :, 34));

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
CairoMakie.activate!()
CairoMakie.save("figures/ke_test_quarter.png", fig, px_per_unit = 5)
=#

fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/16∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:20:70, yticksize = 15, xticksize = 15)
        
hm1 = heatmap!(ax, xF, yF, interior(ζ, :, :, 1)*10^5, colormap = :berlin, colorrange = (-6, 6))
#cb  = Colorbar(fig[0, 2], hm1, vertical = false, label = L"\text{Surface Vertical Vorticity [s}^{-1} \cdot 10^{-5}\text{]}", ticks = ([-3e-5, 0, 3e-5], [L"-3", L"0", L"3"]))
cb = Colorbar(fig[1,2], hm1, width = 30, ticksize = 10, label=L"s^{-1}\cdot 10^{-5}", labelsize = 40, ticklabelsize = 40, height = Relative(3/4))
display(fig)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("figures/vort_sixteen.png", fig, px_per_unit = 5)


fig = Figure(resolution = (1000, 2000))
ax = Axis(fig[1, 1], xlabel="Longitude [∘]", xlabelsize = 40, xticklabelsize = 40, ylabel="Latitude [∘]", ylabelsize = 40,title="1/16∘", yticklabelsize = 40, titlesize=45, aspect=0.5, yticks=-70:20:70, yticksize = 15, xticksize = 15)
hm1 = heatmap!(ax, xC, yC, interior(b, :, :, 34) ./ 2e-3, colormap = :thermal, colorrange = (0, 30))
#cb  = Colorbar(fig[0, 3], hm1, vertical = false, label = L"\text{Surface Temperature [}^\circ\text{C}^{-1}\text{]}", ticks = ([0, 10, 20, 30], [L"0", L"10", L"20", L"30"]))
cb = Colorbar(fig[1,2], hm1, width = 30, height = Relative(3/4), ticksize = 10, label=L"\circ C^{-1}", labelsize = 40, ticklabelsize = 40)
display(fig)
using CairoMakie
CairoMakie.activate!()
CairoMakie.save("figures/surf_temp_sixteen.png", fig, px_per_unit = 5)
